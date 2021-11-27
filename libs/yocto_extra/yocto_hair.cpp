#include "yocto_hair.h"


// -----------------------------------------------------------------------------
// SHADING FUNCTIONS IMPLEMENTATION
// -----------------------------------------------------------------------------

namespace yocto::hair {

hair_data get_hair_data(const material_data& material, float v,
    const vec3f& normal, const vec3f& tangent) {
  
  auto h_data = hair_data{};

  if (material.sigma_a != zero3f) {
    h_data.sigma_a = material.sigma_a;
  } else if (material.color != zero3f) {
    h_data.sigma_a = sigma_a_from_reflectance(material.color, material.beta_n);
  } else if (material.eumelanin || material.pheomelanin) {
    h_data.sigma_a = sigma_a_from_concentration(
        material.eumelanin, material.pheomelanin);
  }

  h_data.eta   = material.eta;
  auto beta_m = material.beta_m;
  auto beta_n = material.beta_n;
  h_data.alpha = material.alpha;
  h_data.h = -1 + 2 * v;

  h_data.gamma_o = safe_asin(h_data.h);

  // Compute longitudinal variance from beta_m
  h_data.v[0] = sqr(
      0.726f * beta_m + 0.812f * sqr(beta_m) + 3.7f * pow<20>(beta_m));
  h_data.v[1] = .25f * h_data.v[0];
  h_data.v[2] = 4 * h_data.v[0];
  for (auto p = 3; p <= p_max; p++) h_data.v[p] = h_data.v[2];

  // Compute azimuthal logistic scale factor from beta_n
  h_data.s = sqrt_pi_over_8f *
           (0.265f * beta_n + 1.194f * sqr(beta_n) + 5.372f * pow<22>(beta_n));

  // Compute alpha terms for hair scales
  h_data.sin_2k_alpha[0] = sin(pif / 180 * h_data.alpha);
  h_data.cos_2k_alpha[0] = safe_sqrt(1 - sqr(h_data.sin_2k_alpha[0]));
  for (auto i = 1; i < 3; i++) {
    h_data.sin_2k_alpha[i] = 2 * h_data.cos_2k_alpha[i - 1] *
                             h_data.sin_2k_alpha[i - 1];
    h_data.cos_2k_alpha[i] = sqr(h_data.cos_2k_alpha[i - 1]) -
                             sqr(h_data.sin_2k_alpha[i - 1]);
  }

  h_data.world_to_bsdf = inverse(frame_fromzy(zero3f, normal, tangent));

  return h_data;
}

vec3f eval_hair_scattering(const hair_data& hair_data, const vec3f& normal, const vec3f& outgoing_,
    const vec3f& incoming_) {
  auto sigma_a       = hair_data.sigma_a;
  auto eta           = hair_data.eta;
  auto h             = hair_data.h;
  auto gamma_o       = hair_data.gamma_o;
  auto v             = hair_data.v;
  auto s             = hair_data.s;
  auto sin_2k_alpha  = hair_data.sin_2k_alpha;
  auto cos_2k_alpha  = hair_data.cos_2k_alpha;
  auto world_to_bsdf = hair_data.world_to_bsdf;

  auto outgoing = transform_direction(world_to_bsdf, outgoing_);
  auto incoming = transform_direction(world_to_bsdf, incoming_);

  // Compute hair coordinate system terms related to outgoing direction
  auto sin_theta_o = outgoing.x;
  auto cos_theta_o = safe_sqrt(1 - sqr(sin_theta_o));
  auto phi_o       = atan2(outgoing.z, outgoing.y);

  // Compute hair coordinate system terms related to incoming direction
  auto sin_theta_i = incoming.x;
  auto cos_theta_i = safe_sqrt(1 - sqr(sin_theta_i));
  auto phi_i       = atan2(incoming.z, incoming.y);

  // Compute cos( theta_t ) for refracted ray
  auto sin_theta_t = sin_theta_o / eta;
  auto cos_theta_t = safe_sqrt(1 - sqr(sin_theta_t));

  // Compute gamma_t for refracted ray
  auto etap        = sqrt(eta * eta - sqr(sin_theta_o)) / cos_theta_o;
  auto sin_gamma_t = h / etap;
  auto cos_gamma_t = safe_sqrt(1 - sqr(sin_gamma_t));
  auto gamma_t     = safe_asin(sin_gamma_t);

  // Compute the transmittance T of a single path through the cylinder
  auto T = exp(-sigma_a * (2 * cos_gamma_t / cos_theta_t));

  // Evaluate hair BSDF
  auto phi  = phi_i - phi_o;
  auto ap   = hair::ap(cos_theta_o, eta, h, T);
  auto fsum = zero3f;

  for (auto p = 0; p < p_max; p++) {
    // Compute sin( theta_o ) and cos( theta_o ) terms accounting for scales
    auto sin_theta_op = 0.0f;
    auto cos_theta_op = 0.0f;
    if (p == 0) {
      sin_theta_op = sin_theta_o * cos_2k_alpha[1] -
                     cos_theta_o * sin_2k_alpha[1];
      cos_theta_op = cos_theta_o * cos_2k_alpha[1] +
                     sin_theta_o * sin_2k_alpha[1];
    }

    // Handle remainder of p values for hair scale tilt
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

    // Handle out-of-range cos( theta_o ) from scale adjustment
    cos_theta_op = abs(cos_theta_op);
    fsum += mp(cos_theta_i, cos_theta_op, sin_theta_i, sin_theta_op, v[p]) *
            ap[p] * np(phi, p, s, gamma_o, gamma_t);
  }

  // Compute contribution of remaining terms after pMax
  fsum += mp(cos_theta_i, cos_theta_o, sin_theta_i, sin_theta_o, v[p_max]) *
          ap[p_max] / (2 * pif);
  return fsum * abs(dot(normal, incoming));
}

vec3f sample_hair_scattering(
    const hair_data& hair_data, const vec3f& outgoing_, const vec2f& rn) {
  auto eta           = hair_data.eta;
  auto h             = hair_data.h;
  auto gamma_o       = hair_data.gamma_o;
  auto v             = hair_data.v;
  auto s             = hair_data.s;
  auto sin_2k_alpha  = hair_data.sin_2k_alpha;
  auto cos_2k_alpha  = hair_data.cos_2k_alpha;

  auto outgoing = transform_direction(hair_data.world_to_bsdf, outgoing_);

  // Compute hair coordinate system terms related to outgoing direction
  auto sin_theta_o = outgoing.x;
  auto cos_theta_o = safe_sqrt(1 - sqr(sin_theta_o));
  auto phi_o       = atan2(outgoing.z, outgoing.y);

  // Derive four random samples from rn
  auto u = std::array<vec2f, 2>{demux_float(rn.x), demux_float(rn.y)};

  // Determine which term p to sample for hair scattering
  auto ap_pdf = compute_ap_pdf(hair_data, cos_theta_o);
  auto p      = 0;
  for (p = 0; p < p_max; p++) {
    if (u[0][0] < ap_pdf[p]) break;
    u[0][0] -= ap_pdf[p];
  }

  // Rotate sin( theta_o ) and cos( theta_o ) to account for hair scale tilt
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

  // Sample M_p to compute theta_i
  u[1][0]          = max(u[1][0], 1e-5f);
  auto cos_theta   = 1 + v[p] * log(u[1][0] + (1 - u[1][0]) * exp(-2 / v[p]));
  auto sin_theta   = safe_sqrt(1 - sqr(cos_theta));
  auto cos_phi     = cos(2 * pif * u[1][1]);
  auto sin_theta_i = -cos_theta * sin_theta_op +
                     sin_theta * cos_phi * cos_theta_op;
  auto cos_theta_i = safe_sqrt(1 - sqr(sin_theta_i));

  // Sample N_p to compute delta phi

  // Compute gamma_t for refracted ray
  auto etap        = sqrt(eta * eta - sqr(sin_theta_o)) / cos_theta_o;
  auto sin_gamma_t = h / etap;
  auto gamma_t     = safe_asin(sin_gamma_t);
  auto dphi        = 0.0f;
  if (p < p_max)
    dphi = phi(p, gamma_o, gamma_t) +
           sample_trimmed_logistic(u[0][1], s, -pif, pif);
  else
    dphi = 2 * pif * u[0][1];

  // Compute incoming direction from sampled hair scattering angles
  auto phi_i = phi_o + dphi;

  auto incoming = vec3f{
      sin_theta_i, cos_theta_i * cos(phi_i), cos_theta_i * sin(phi_i)};
  
  return transform_direction(inverse(hair_data.world_to_bsdf), incoming);
}


float sample_hair_scattering_pdf(const hair_data& hair_data,
    const vec3f& outgoing_, const vec3f& incoming_) {
  auto eta           = hair_data.eta;
  auto h             = hair_data.h;
  auto gamma_o       = hair_data.gamma_o;
  auto v             = hair_data.v;
  auto s             = hair_data.s;
  auto sin_2k_alpha  = hair_data.sin_2k_alpha;
  auto cos_2k_alpha  = hair_data.cos_2k_alpha;

  auto outgoing = transform_direction(hair_data.world_to_bsdf, outgoing_);
  auto incoming = transform_direction(hair_data.world_to_bsdf, incoming_);

  // Compute hair coordinate system terms related to outgoing direction
  auto sin_theta_o = outgoing.x;
  auto cos_theta_o = safe_sqrt(1 - sqr(sin_theta_o));
  auto phi_o       = atan2(outgoing.z, outgoing.y);

  // Compute hair coordinate system terms related to incoming direction
  auto sin_theta_i = incoming.x;
  auto cos_theta_i = safe_sqrt(1 - sqr(sin_theta_i));
  auto phi_i       = atan2(incoming.z, incoming.y);

  // Compute gamma_t for refracted ray
  auto etap        = sqrt(eta * eta - sqr(sin_theta_o)) / cos_theta_o;
  auto sin_gamma_t = h / etap;
  auto gamma_t     = safe_asin(sin_gamma_t);

  // Compute PDF for A_p terms
  auto ap_pdf = compute_ap_pdf(hair_data, cos_theta_o);

  // Compute PDF sum for hair scattering events
  auto phi = phi_i - phi_o;
  auto pdf = 0.0f;
  for (auto p = 0; p < p_max; p++) {
    // Compute sin( theta_o ) and cos( theta_o ) terms accounting for scales
    auto sin_theta_op = 0.0f;
    auto cos_theta_op = 0.0f;
    if (p == 0) {
      sin_theta_op = sin_theta_o * cos_2k_alpha[1] -
                     cos_theta_o * sin_2k_alpha[1];
      cos_theta_op = cos_theta_o * cos_2k_alpha[1] +
                     sin_theta_o * sin_2k_alpha[1];
    }

    // Handle remainder of p values for hair scale tilt
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

    // Handle out-of-range cos( theta_o ) from scale adjustment
    cos_theta_op = abs(cos_theta_op);
    pdf += mp(cos_theta_i, cos_theta_op, sin_theta_i, sin_theta_op, v[p]) *
           ap_pdf[p] * np(phi, p, s, gamma_o, gamma_t);
  }
  pdf += mp(cos_theta_i, cos_theta_o, sin_theta_i, sin_theta_o, v[p_max]) *
         ap_pdf[p_max] * (1 / (2 * pif));

  return pdf;
}

}  // namespace yocto::hair