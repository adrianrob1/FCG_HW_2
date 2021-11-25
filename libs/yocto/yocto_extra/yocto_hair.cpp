#include "yocto_hair.h"


// -----------------------------------------------------------------------------
// SHADING FUNCTIONS IMPLEMENTATION
// -----------------------------------------------------------------------------

namespace yocto::hair {

hair_brdf eval_hair_brdf(const hair_material& material, float v,
    const vec3f& normal, const vec3f& tangent) {
  
  auto brdf = hair_brdf{};

  if (material.sigma_a != zero3f) {
    brdf.sigma_a = material.sigma_a;
  } else if (material.color != zero3f) {
    brdf.sigma_a = sigma_a_from_reflectance(material.color, material.beta_n);
  } else if (material.eumelanin || material.pheomelanin) {
    brdf.sigma_a = sigma_a_from_concentration(
        material.eumelanin, material.pheomelanin);
  }

  auto beta_m = material.beta_m;
  auto beta_n = material.beta_n;
  brdf.alpha  = material.alpha;
  brdf.eta    = material.eta;

  brdf.h = -1 + 2 * v;

  brdf.gamma_o = safe_asin(brdf.h);

  // Compute longitudinal variance from beta_m
  brdf.v[0] = sqr(
      0.726f * beta_m + 0.812f * sqr(beta_m) + 3.7f * pow<20>(beta_m));
  brdf.v[1] = 0.25f * brdf.v[0];
  brdf.v[2] = 4 * brdf.v[0];
  for (auto p = 3; p <= p_max; p++) brdf.v[p] = brdf.v[2];

  // Compute azimuthal logistic scale factor from beta_n
  brdf.s = sqrt_pi_over_8f *
           (0.265f * beta_n + 1.194f * sqr(beta_n) + 5.372f * pow<22>(beta_n));

  // Compute alpha terms for hair scales
  brdf.sin_2k_alpha[0] = sin(pif / 180 * brdf.alpha);
  brdf.cos_2k_alpha[0] = safe_sqrt(1 - sqr(brdf.sin_2k_alpha[0]));
  for (auto i = 1; i < 3; i++) {
    brdf.sin_2k_alpha[i] = 2 * brdf.cos_2k_alpha[i - 1] *
                           brdf.sin_2k_alpha[i - 1];
    brdf.cos_2k_alpha[i] = sqr(brdf.cos_2k_alpha[i - 1]) -
                           sqr(brdf.sin_2k_alpha[i - 1]);
  }

  brdf.world_to_brdf = inverse(frame_fromzx(zero3f, normal, tangent));

  return brdf;
}

vec3f eval_hair_scattering(
    const hair_brdf& brdf, const vec3f& outgoing_, const vec3f& incoming_) {
  auto sigma_a       = brdf.sigma_a;
  auto eta           = brdf.eta;
  auto h             = brdf.h;
  auto gamma_o       = brdf.gamma_o;
  auto v             = brdf.v;
  auto s             = brdf.s;
  auto sin_2k_alpha  = brdf.sin_2k_alpha;
  auto cos_2k_alpha  = brdf.cos_2k_alpha;
  auto world_to_brdf = brdf.world_to_brdf;

  auto outgoing = transform_direction(world_to_brdf, outgoing_);
  auto incoming = transform_direction(world_to_brdf, incoming_);

  // Compute hair coordinate system terms related to _wo_
  auto sin_theta_o = outgoing.x;
  auto cos_theta_o = safe_sqrt(1 - sqr(sin_theta_o));
  auto phi_o       = atan2(outgoing.z, outgoing.y);

  // Compute hair coordinate system terms related to _wi_
  auto sin_theta_i = incoming.x;
  auto cos_theta_i = safe_sqrt(1 - sqr(sin_theta_i));
  auto phi_i       = atan2(incoming.z, incoming.y);

  // Compute $\cos \thetat$ for refracted ray
  auto sin_theta_t = sin_theta_o / eta;
  auto cos_theta_t = safe_sqrt(1 - sqr(sin_theta_t));

  // Compute $\gammat$ for refracted ray
  auto etap        = sqrt(eta * eta - sqr(sin_theta_o)) / cos_theta_o;
  auto sin_gamma_t = h / etap;
  auto cos_gamma_t = safe_sqrt(1 - sqr(sin_gamma_t));
  auto gamma_t     = safe_asin(sin_gamma_t);

  // Compute the transmittance _T_ of a single path through the cylinder
  auto T = exp(-sigma_a * (2 * cos_gamma_t / cos_theta_t));

  // Evaluate hair BSDF
  auto phi  = phi_i - phi_o;
  auto ap   = hair::ap(cos_theta_o, eta, h, T);
  auto fsum = zero3f;

  for (auto p = 0; p < p_max; p++) {
    // Compute $\sin \thetao$ and $\cos \thetao$ terms accounting for scales
    auto sin_theta_op = 0.0f;
    auto cos_theta_op = 0.0f;
    if (p == 0) {
      sin_theta_op = sin_theta_o * cos_2k_alpha[1] -
                     cos_theta_o * sin_2k_alpha[1];
      cos_theta_op = cos_theta_o * cos_2k_alpha[1] +
                     sin_theta_o * sin_2k_alpha[1];
    }

    // Handle remainder of $p$ values for hair scale tilt
    else if (p == 1) {
      sin_theta_op = sin_theta_o * cos_2k_alpha[0] +
                     cos_theta_o * sin_2k_alpha[0];
      cos_theta_op = cos_theta_o * cos_2k_alpha[0] -
                     sin_theta_o * sin_2k_alpha[0];
    } else if (p == 2) {
      sin_theta_op = sin_theta_o * cos_2k_alpha[2] +
                     cos_theta_o * sin_2k_alpha[2];
      cos_theta_op = cos_theta_o * cos_2k_alpha[2] -
                     sin_theta_o * sin_2k_alpha[2];
    } else {
      sin_theta_op = sin_theta_o;
      cos_theta_op = cos_theta_o;
    }

    // Handle out-of-range $\cos \thetao$ from scale adjustment
    cos_theta_op = abs(cos_theta_op);
    fsum += mp(cos_theta_i, cos_theta_op, sin_theta_i, sin_theta_op, v[p]) *
            ap[p] * np(phi, p, s, gamma_o, gamma_t);
  }

  // Compute contribution of remaining terms after _pMax_
  fsum += mp(cos_theta_i, cos_theta_o, sin_theta_i, sin_theta_o, v[p_max]) *
          ap[p_max] / (2 * pif);
  // if (abs(incoming.z) > 0)
  //   fsum /= abs(incoming.z);
  return fsum;
}

vec3f sample_hair_scattering(
    const hair_brdf& brdf, const vec3f& outgoing_, const vec2f& rn) {
  auto eta           = brdf.eta;
  auto h             = brdf.h;
  auto gamma_o       = brdf.gamma_o;
  auto v             = brdf.v;
  auto s             = brdf.s;
  auto sin_2k_alpha  = brdf.sin_2k_alpha;
  auto cos_2k_alpha  = brdf.cos_2k_alpha;
  auto world_to_brdf = brdf.world_to_brdf;

  auto outgoing = transform_direction(world_to_brdf, outgoing_);

  // Compute hair coordinate system terms related to _wo_
  auto sin_theta_o = outgoing.x;
  auto cos_theta_o = safe_sqrt(1 - sqr(sin_theta_o));
  auto phi_o       = atan2(outgoing.z, outgoing.y);

  // Derive four random samples from _u2_
  auto u = std::array<vec2f, 2>{demux_float(rn.x), demux_float(rn.y)};

  // Determine which term $p$ to sample for hair scattering
  auto ap_pdf = compute_ap_pdf(brdf, cos_theta_o);
  auto p      = 0;
  for (p = 0; p < p_max; p++) {
    if (u[0][0] < ap_pdf[p]) break;
    u[0][0] -= ap_pdf[p];
  }

  // Rotate $\sin \thetao$ and $\cos \thetao$ to account for hair scale tilt
  auto sin_theta_op = 0.0f;
  auto cos_theta_op = 0.0f;
  if (p == 0) {
    sin_theta_op = sin_theta_o * cos_2k_alpha[1] -
                   cos_theta_o * sin_2k_alpha[1];
    cos_theta_op = cos_theta_o * cos_2k_alpha[1] +
                   sin_theta_o * sin_2k_alpha[1];
  } else if (p == 1) {
    sin_theta_op = sin_theta_o * cos_2k_alpha[0] +
                   cos_theta_o * sin_2k_alpha[0];
    cos_theta_op = cos_theta_o * cos_2k_alpha[0] -
                   sin_theta_o * sin_2k_alpha[0];
  } else if (p == 2) {
    sin_theta_op = sin_theta_o * cos_2k_alpha[2] +
                   cos_theta_o * sin_2k_alpha[2];
    cos_theta_op = cos_theta_o * cos_2k_alpha[2] -
                   sin_theta_o * sin_2k_alpha[2];
  } else {
    sin_theta_op = sin_theta_o;
    cos_theta_op = cos_theta_o;
  }

  // Sample $M_p$ to compute $\thetai$
  u[1][0]          = max(u[1][0], 1e-5f);
  auto cos_theta   = 1 + v[p] * log(u[1][0] + (1 - u[1][0]) * exp(-2 / v[p]));
  auto sin_theta   = safe_sqrt(1 - sqr(cos_theta));
  auto cos_phi     = cos(2 * pif * u[1][1]);
  auto sin_theta_i = -cos_theta * sin_theta_op +
                     sin_theta * cos_phi * cos_theta_op;
  auto cos_theta_i = safe_sqrt(1 - sqr(sin_theta_i));

  // Sample $N_p$ to compute $\Delta\phi$

  // Compute $\gammat$ for refracted ray
  auto etap        = sqrt(eta * eta - sqr(sin_theta_o)) / cos_theta_o;
  auto sin_gamma_t = h / etap;
  auto gamma_t     = safe_asin(sin_gamma_t);
  auto dphi        = 0.0f;
  if (p < p_max)
    dphi = phi(p, gamma_o, gamma_t) +
           sample_trimmed_logistic(u[0][1], s, -pif, pif);
  else
    dphi = 2 * pif * u[0][1];

  // Compute _wi_ from sampled hair scattering angles
  auto phi_i = phi_o + dphi;

  auto incoming = vec3f{
      sin_theta_i, cos_theta_i * cos(phi_i), cos_theta_i * sin(phi_i)};
  return transform_direction(inverse(world_to_brdf), incoming);
}
float sample_hair_scattering_pdf(
    const hair_brdf& brdf, const vec3f& outgoing_, const vec3f& incoming_) {
  auto eta           = brdf.eta;
  auto h             = brdf.h;
  auto gamma_o       = brdf.gamma_o;
  auto v             = brdf.v;
  auto s             = brdf.s;
  auto sin_2k_alpha  = brdf.sin_2k_alpha;
  auto cos_2k_alpha  = brdf.cos_2k_alpha;
  auto world_to_brdf = brdf.world_to_brdf;

  auto outgoing = transform_direction(world_to_brdf, outgoing_);
  auto incoming = transform_direction(world_to_brdf, incoming_);

  // Compute hair coordinate system terms related to _wo_
  auto sin_theta_o = outgoing.x;
  auto cos_theta_o = safe_sqrt(1 - sqr(sin_theta_o));
  auto phi_o       = atan2(outgoing.z, outgoing.y);

  // Compute hair coordinate system terms related to _wi_
  auto sin_theta_i = incoming.x;
  auto cos_theta_i = safe_sqrt(1 - sqr(sin_theta_i));
  auto phi_i       = atan2(incoming.z, incoming.y);

  // Compute $\gammat$ for refracted ray
  auto etap        = sqrt(eta * eta - sqr(sin_theta_o)) / cos_theta_o;
  auto sin_gamma_t = h / etap;
  auto gamma_t     = safe_asin(sin_gamma_t);

  // Compute PDF for $A_p$ terms
  auto ap_pdf = compute_ap_pdf(brdf, cos_theta_o);

  // Compute PDF sum for hair scattering events
  auto phi = phi_i - phi_o;
  auto pdf = 0.0f;
  for (auto p = 0; p < p_max; p++) {
    // Compute $\sin \thetao$ and $\cos \thetao$ terms accounting for scales
    auto sin_theta_op = 0.0f;
    auto cos_theta_op = 0.0f;
    if (p == 0) {
      sin_theta_op = sin_theta_o * cos_2k_alpha[1] -
                     cos_theta_o * sin_2k_alpha[1];
      cos_theta_op = cos_theta_o * cos_2k_alpha[1] +
                     sin_theta_o * sin_2k_alpha[1];
    }

    // Handle remainder of $p$ values for hair scale tilt
    else if (p == 1) {
      sin_theta_op = sin_theta_o * cos_2k_alpha[0] +
                     cos_theta_o * sin_2k_alpha[0];
      cos_theta_op = cos_theta_o * cos_2k_alpha[0] -
                     sin_theta_o * sin_2k_alpha[0];
    } else if (p == 2) {
      sin_theta_op = sin_theta_o * cos_2k_alpha[2] +
                     cos_theta_o * sin_2k_alpha[2];
      cos_theta_op = cos_theta_o * cos_2k_alpha[2] -
                     sin_theta_o * sin_2k_alpha[2];
    } else {
      sin_theta_op = sin_theta_o;
      cos_theta_op = cos_theta_o;
    }

    // Handle out-of-range $\cos \thetao$ from scale adjustment
    cos_theta_op = abs(cos_theta_op);
    pdf += mp(cos_theta_i, cos_theta_op, sin_theta_i, sin_theta_op, v[p]) *
           ap_pdf[p] * np(phi, p, s, gamma_o, gamma_t);
  }
  pdf += mp(cos_theta_i, cos_theta_o, sin_theta_i, sin_theta_o, v[p_max]) *
         ap_pdf[p_max] * (1 / (2 * pif));
  return pdf;
}

}  // namespace yocto::hair