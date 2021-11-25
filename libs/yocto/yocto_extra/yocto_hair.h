#pragma once

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <../yocto_math.h>
#include <../yocto_shading.h>

// -----------------------------------------------------------------------------
// MATERIALS
// -----------------------------------------------------------------------------

namespace yocto::hair {

inline const int p_max = 3;

struct hair_material {
  vec3f sigma_a     = zero3f;   // absorption coefficient of the interior
  float beta_m      = 0.3;      // longitudinal roughness (along length)
  float beta_n      = 0.3;      // azimuthal roughness (along width)
  float alpha       = 2;        // angle of hair's scales
  float eta         = 1.55;     // ior of the interior
  vec3f color       = zero3f;   // hair color
  float eumelanin   = 0;
  float pheomelanin = 0;
};

struct hair_brdf {
  vec3f sigma_a = zero3f;
  float alpha   = 2;
  float eta     = 1.55;
  float h       = 0;

  // computed properties
  std::array<float, p_max + 1> v;
  float                        s = 0;
  vec3f                        sin_2k_alpha;
  vec3f                        cos_2k_alpha;
  float                        gamma_o = 0;

  // Convert outgoing and incoming directions to BRDF coordinate system
  frame3f world_to_brdf = identity3x4f;
};

}  // namespace yocto::hair

// -----------------------------------------------------------------------------
// SHADING FUNCTIONS
// -----------------------------------------------------------------------------

namespace yocto::hair {

hair_brdf eval_hair_brdf(const hair_material& material, float v,
    const vec3f& normal, const vec3f& tangent);

vec3f eval_hair_scattering(
    const hair_brdf& brdf, const vec3f& outgoing, const vec3f& incoming);

vec3f sample_hair_scattering(
    const hair_brdf& brdf, const vec3f& outgoing, const vec2f& rn);

float sample_hair_scattering_pdf(
    const hair_brdf& brdf, const vec3f& outgoing, const vec3f& incoming);

}  // namespace yocto::hair


// -----------------------------------------------------------------------------
// MATH WORKER FUNCTIONS
// -----------------------------------------------------------------------------

namespace yocto::hair {

inline const float sqrt_pi_over_8f = 0.626657069f;

template <int n>
static float pow(float v) {
  auto n2 = pow<n / 2>(v);
  return n2 * n2 * pow<n & 1>(v);
}

template <>
inline float pow<1>(float v) {
  return v;
}

template <>
inline float pow<0>(float v) {
  return 1;
}

inline float safe_asin(float x) { return asin(clamp(x, -1.0f, 1.0f)); }

inline float safe_sqrt(float x) { return sqrt(max(0.0f, x)); }

static vec3f sigma_a_from_concentration(float ce, float cp) {
  auto eumelanin_sigma_a   = vec3f{0.419f, 0.697f, 1.37f};
  auto pheomelanin_sigma_a = vec3f{0.187f, 0.4f, 1.05f};
  return ce * eumelanin_sigma_a + cp * pheomelanin_sigma_a;
}

static vec3f sigma_a_from_reflectance(const vec3f& c, float beta_n) {
  return sqr(log(c) / (5.969f - 0.215f * beta_n + 2.532f * sqr(beta_n) -
                          10.73f * pow<3>(beta_n) + 5.574f * pow<4>(beta_n) +
                          0.245f * pow<5>(beta_n)));
}

inline float i0(float x) {
  float   val   = 0;
  float   x2i   = 1;
  int64_t ifact = 1;
  int     i4    = 1;
  // I0(x) \approx Sum_i x^(2i) / (4^i (i!)^2)
  for (int i = 0; i < 10; i++) {
    if (i > 1) ifact *= i;
    val += x2i / (i4 * ifact * ifact);
    x2i *= x * x;
    i4 *= 4;
  }
  return val;
}

inline float log_i0(float x) {
  if (x > 12)
    return x + 0.5f * (-log(2 * pif) + log(1 / x) + 1 / (8 * x));
  else
    return log(i0(x));
}

static float mp(float cos_theta_i, float cos_theta_o, float sin_theta_i,
    float sin_theta_o, float v) {
  auto a = cos_theta_i * cos_theta_o / v;
  auto b = sin_theta_i * sin_theta_o / v;
  return (v <= 0.1f) ? (exp(log_i0(a) - b - 1 / v + 0.6931f + log(1 / (2 * v))))
                     : (exp(-b) * i0(a)) / (sinh(1 / v) * 2 * v);
}

static std::array<vec3f, p_max + 1> ap(
    float cos_theta_o, float eta, float h, const vec3f& T) {
  auto ap = std::array<vec3f, p_max + 1>{};
  // Compute $p=0$ attenuation at initial cylinder intersection
  auto cos_gamma_o = safe_sqrt(1 - h * h);
  auto cos_theta   = cos_theta_o * cos_gamma_o;

  // fresnel_dielectric(eta, cos_theta): we hack function's parameters bulding
  // two vectors s.t. their dot product gives exactly cos_theta
  auto f = fresnel_dielectric(eta, {0, 0, 1}, {0, 0, cos_theta});
  ap[0]  = {f, f, f};

  // Compute $p=1$ attenuation term
  ap[1] = sqr(1 - f) * T;

  // Compute attenuation terms up to $p=_pMax_$
  for (auto p = 2; p < p_max; p++) ap[p] = ap[p - 1] * T * f;

  // Compute attenuation term accounting for remaining orders of scattering
  ap[p_max] = ap[p_max - 1] * f * T / (vec3f{1, 1, 1} - T * f);
  return ap;
}

inline float phi(int p, float gamma_o, float gamma_t) {
  return 2 * p * gamma_t - 2 * gamma_o + p * pif;
}

inline float logistic(float x, float s) {
  x = abs(x);
  return exp(-x / s) / (s * sqr(1 + exp(-x / s)));
}

inline float logistic_cdf(float x, float s) { return 1 / (1 + exp(-x / s)); }

inline float trimmed_logistic(float x, float s, float a, float b) {
  return logistic(x, s) / (logistic_cdf(b, s) - logistic_cdf(a, s));
}

inline float np(float phi, int p, float s, float gamma_o, float gamma_t) {
  auto dphi = phi - hair::phi(p, gamma_o, gamma_t);
  // Remap _dphi_ to $[-\pi,\pi]$
  while (dphi > pif) dphi -= 2 * pif;
  while (dphi < -pif) dphi += 2 * pif;
  return trimmed_logistic(dphi, s, -pif, pif);
}

// https://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/
static uint32_t compact1by1(uint32_t x) {
  // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
  x &= 0x55555555;
  // x = --fe --dc --ba --98 --76 --54 --32 --10
  x = (x ^ (x >> 1)) & 0x33333333;
  // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
  x = (x ^ (x >> 2)) & 0x0f0f0f0f;
  // x = ---- ---- fedc ba98 ---- ---- 7654 3210
  x = (x ^ (x >> 4)) & 0x00ff00ff;
  // x = ---- ---- ---- ---- fedc ba98 7654 3210
  x = (x ^ (x >> 8)) & 0x0000ffff;
  return x;
}

static vec2f demux_float(float f) {
  uint64_t v       = f * (1ull << 32);
  uint32_t bits[2] = {compact1by1(v), compact1by1(v >> 1)};
  return {bits[0] / float(1 << 16), bits[1] / float(1 << 16)};
}

static float sample_trimmed_logistic(float u, float s, float a, float b) {
  auto k = logistic_cdf(b, s) - logistic_cdf(a, s);
  auto x = -s * log(1 / (u * k + logistic_cdf(a, s)) - 1);
  return clamp(x, a, b);
}

// Luminance
inline float luminance(const vec3f& a) {
  return (0.2126f * a.x + 0.7152f * a.y + 0.0722f * a.z);
}

static std::array<float, p_max + 1> compute_ap_pdf(
    const hair_brdf& brdf, float cos_theta_o) {
  auto sigma_a = brdf.sigma_a;
  auto eta     = brdf.eta;
  auto h       = brdf.h;

  // Compute array of $A_p$ values for _cosThetaO_
  auto sin_theta_o = safe_sqrt(1 - cos_theta_o * cos_theta_o);

  // Compute $\cos \thetat$ for refracted ray
  auto sin_theta_t = sin_theta_o / eta;
  auto cos_theta_t = safe_sqrt(1 - sqr(sin_theta_t));

  // Compute $\gammat$ for refracted ray
  auto etap        = sqrt(eta * eta - sqr(sin_theta_o)) / cos_theta_o;
  auto sin_gamma_t = h / etap;
  auto cos_gamma_t = safe_sqrt(1 - sqr(sin_gamma_t));

  // Compute the transmittance _T_ of a single path through the cylinder
  auto T  = exp(-sigma_a * (2 * cos_gamma_t / cos_theta_t));
  auto ap = hair::ap(cos_theta_o, eta, h, T);

  // Compute $A_p$ PDF from individual $A_p$ terms
  auto ap_pdf = std::array<float, p_max + 1>{};
  auto sum_y  = 0.0f;
  for (auto i = 0; i <= p_max; i++) {
    sum_y += luminance(ap[i]);
  }
  for (auto i = 0; i <= p_max; i++) {
    ap_pdf[i] = luminance(ap[i]) / sum_y;
  }
  return ap_pdf;
}

}