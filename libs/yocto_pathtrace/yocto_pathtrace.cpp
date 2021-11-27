//
// Implementation for Yocto/PathTrace.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

#include "yocto_pathtrace.h"

#include <yocto/yocto_cli.h>
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_parallel.h>
#include <yocto/yocto_sampling.h>
#include <yocto/yocto_shading.h>
#include <yocto/yocto_shape.h>

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PATH TRACING
// -----------------------------------------------------------------------------
namespace yocto {

// Convenience functions
[[maybe_unused]] static vec3f eval_position(
    const scene_data& scene, const bvh_intersection& intersection) {
  return eval_position(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv);
}
[[maybe_unused]] static vec3f eval_normal(
    const scene_data& scene, const bvh_intersection& intersection) {
  return eval_normal(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv);
}
[[maybe_unused]] static vec3f eval_element_normal(
    const scene_data& scene, const bvh_intersection& intersection) {
  return eval_element_normal(
      scene, scene.instances[intersection.instance], intersection.element);
}
[[maybe_unused]] static vec3f eval_shading_position(const scene_data& scene,
    const bvh_intersection& intersection, const vec3f& outgoing) {
  return eval_shading_position(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv, outgoing);
}
[[maybe_unused]] static vec3f eval_shading_normal(const scene_data& scene,
    const bvh_intersection& intersection, const vec3f& outgoing) {
  return eval_shading_normal(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv, outgoing);
}
[[maybe_unused]] static vec2f eval_texcoord(
    const scene_data& scene, const bvh_intersection& intersection) {
  return eval_texcoord(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv);
}
[[maybe_unused]] static material_point eval_material(
    const scene_data& scene, const bvh_intersection& intersection) {
  return eval_material(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv);
}
[[maybe_unused]] static bool is_volumetric(
    const scene_data& scene, const bvh_intersection& intersection) {
  return is_volumetric(scene, scene.instances[intersection.instance]);
}

// Evaluates/sample the BRDF scaled by the cosine of the incoming direction.
static vec3f eval_emission(const material_point& material, const vec3f& normal,
    const vec3f& outgoing) {
  return dot(normal, outgoing) >= 0 ? material.emission : vec3f{0, 0, 0};
}

// Evaluates/sample the BRDF scaled by the cosine of the incoming direction.
static vec3f eval_bsdfcos(const material_point& material, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
  // YOUR CODE GOES HERE

  if (material.roughness == 0 && material.type != material_type::hair)
    return {0, 0, 0};

  switch (material.type) {
    case material_type::matte:
      return eval_matte(material.color, normal, outgoing, incoming);
    case material_type::glossy:
      return eval_glossy(material.color, material.ior, material.roughness,
          normal, outgoing, incoming);
    case material_type::reflective:
      return eval_reflective(
          material.color, material.roughness, normal, outgoing, incoming);
    case material_type::transparent:
      return eval_transparent(material.color, material.ior, material.roughness,
          normal, outgoing, incoming);
    case material_type::refractive:
      return eval_refractive(material.color, material.ior, material.roughness,
          normal, outgoing, incoming);
    case material_type::hair:
      return hair::eval_hair_scattering(material.hair, outgoing, incoming);
    default: return {0, 0, 0};
  }
}

static vec3f eval_delta(const material_point& material, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
  // YOUR CODE GOES HERE
  if (material.roughness != 0) return {0, 0, 0};

  switch (material.type) {
    case material_type::reflective:
      return eval_reflective(material.color, normal, outgoing, incoming);
    case material_type::transparent:
      return eval_transparent(
          material.color, material.ior, normal, outgoing, incoming);
    case material_type::refractive:
      return eval_refractive(
          material.color, material.ior, normal, outgoing, incoming);
    default: return {0, 0, 0};
  }
}

// Picks a direction based on the BRDF
static vec3f sample_bsdfcos(const material_point& material, const vec3f& normal,
    const vec3f& outgoing, float rnl, const vec2f& rn) {
  // YOUR CODE GOES HERE
  if (material.roughness == 0 && material.type != material_type::hair)
    return {0, 0, 0};

  switch (material.type) {
    case material_type::matte:
      return sample_matte(material.color, normal, outgoing, rn);
    case material_type::glossy:
      return sample_glossy(material.color, material.ior, material.roughness,
          normal, outgoing, rnl, rn);
    case material_type::reflective:
      return sample_reflective(
          material.color, material.roughness, normal, outgoing, rn);
    case material_type::transparent:
      return sample_transparent(material.color, material.ior,
          material.roughness, normal, outgoing, rnl, rn);
    case material_type::refractive:
      return sample_refractive(material.color, material.ior, material.roughness,
          normal, outgoing, rnl, rn);
    case material_type::hair:
      return hair::sample_hair_scattering(material.hair, outgoing, rn);
    default: return {0, 0, 0};
  }
}

static vec3f sample_delta(const material_point& material, const vec3f& normal,
    const vec3f& outgoing, float rnl) {
  // YOUR CODE GOES HERE

  if (material.roughness != 0) return {0, 0, 0};

  switch (material.type) {
    case material_type::reflective:
      return sample_reflective(material.color, normal, outgoing);
    case material_type::transparent:
      return sample_transparent(
          material.color, material.ior, normal, outgoing, rnl);
    case material_type::refractive:
      return sample_refractive(
          material.color, material.ior, normal, outgoing, rnl);
    default: return {0, 0, 0};
  }
}

// Compute the weight for sampling the BRDF
static float sample_bsdfcos_pdf(const material_point& material,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
  // YOUR CODE GOES HERE

   if (material.roughness == 0 && material.type != material_type::hair)
    return 0;

  switch (material.type) {
    case material_type::matte:
      return sample_matte_pdf(material.color, normal, outgoing, incoming);
    case material_type::glossy:
      return sample_glossy_pdf(material.color, material.ior, material.roughness,
          normal, outgoing, incoming);
    case material_type::reflective:
      return sample_reflective_pdf(
          material.color, material.roughness, normal, outgoing, incoming);
    case material_type::transparent:
      return sample_tranparent_pdf(material.color, material.ior,
          material.roughness, normal, outgoing, incoming);
    case material_type::refractive:
      return sample_refractive_pdf(material.color, material.ior,
          material.roughness, normal, outgoing, incoming);
    case material_type::hair:
      return hair::sample_hair_scattering_pdf(material.hair, outgoing, incoming);
    default: return 0;
  }
}

static float sample_delta_pdf(const material_point& material,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
  // YOUR CODE GOES HERE

  if (material.roughness != 0) return 0;

  switch (material.type) {
    case material_type::reflective:
      return sample_reflective_pdf(material.color, normal, outgoing, incoming);
    case material_type::transparent:
      return sample_tranparent_pdf(
          material.color, material.ior, normal, outgoing, incoming);
    case material_type::refractive:
      return sample_refractive_pdf(
          material.color, material.ior, normal, outgoing, incoming);
    default: return 0;
  }
}

// sample using alias method
inline int sample_discrete_alias(
    const vector<double>& probs, const vector<int>& alias, float ru, float rb) {
  // pick a column with uniform probability
  auto idx = sample_uniform(probs.size(), ru);

  // toss biased coin to determine which option to pick
  auto coinToss = rb < probs[idx];

  return coinToss ? idx : alias[idx];
}

// Sample lights wrt solid angle
static vec3f sample_lights(const scene_data& scene,
    const pathtrace_lights& lights, const vec3f& position, float rl, float rel,
    float rb, const vec2f& ruv) {
  // YOUR CODE GOES HERE
  auto  lid   = sample_uniform((int)lights.lights.size(), rl);
  auto& light = lights.lights[lid];

  if (light.instance != invalidid) {
    auto& instance = scene.instances[light.instance];
    auto& shape    = scene.shapes[instance.shape];
    auto  uv       = ruv;
    if (!shape.triangles.empty()) uv = sample_triangle(ruv);

    // TODO: change to next random triangle id
    /*auto tid = sample_discrete(light.elements_cdf, rel);*/
    auto tid = sample_discrete_alias(
        light.probabilities_modified, light.alias, rel, rb);
    auto lp = eval_position(scene, instance, tid, uv);

    return normalize(lp - position);
  } else if (light.environment != invalidid) {
    auto& env     = scene.environments[light.environment];
    auto& texture = scene.textures[env.emission_tex];

    // TODO: change to next random texel id
    /*auto tid = sample_discrete(light.elements_cdf, rel);*/
    auto tid = sample_discrete_alias(
        light.probabilities_modified, light.alias, rel, rb);

    auto uv = vec2f{((tid % texture.width) + 0.5f) / texture.width,
        ((tid / texture.width) + 0.5f) / texture.height};

    return transform_direction(
        env.frame, {cosf(uv.x * 2 * pif) * sinf(uv.y * pif), cosf(uv.y * pif),
                       sinf(uv.x * 2 * pif) * sinf(uv.y * pif)});
  }

  return {0, 0, 0};
}

// Sample lights pdf
static float sample_lights_pdf(const scene_data& scene, const bvh_data& bvh,
    const pathtrace_lights& lights, const vec3f& position,
    const vec3f& direction) {
  // YOUR CODE GOES HERE
  // the pdf is the sum of all the possible ways
  // that a light ray gets generated

  auto pdf      = 0.0f;
  auto next_pos = position;
  for (auto& light : lights.lights) {
    if (light.instance != invalidid) {
      auto& instance = scene.instances[light.instance];
      auto  lpdf     = 0.0f;

      for (auto bounce = 0; bounce < 100; bounce++) {
        // sum all the possibilities that I pick that direction
        // sum over all the possible intersections from the shading position
        // and accumulate the probability of that one
        auto isec = intersect_bvh(
            bvh, scene, light.instance, {next_pos, direction});
        if (!isec.hit) break;

        // compute pdf on this triangle
        auto lposition = eval_position(scene, instance, isec.element, isec.uv);
        //auto lnormal   = eval_element_normal(scene, instance, isec.element);
        auto lnormal   = eval_normal(scene, instance, isec.element, isec.uv);

        // the last element contains the area of the light
        /*auto area = light.elements_cdf.back();*/
        auto area = light.area;

        // accumulate the probability
        lpdf += distance_squared(lposition, position) /
                (abs(dot(lnormal, direction)) * area);

        // continue
        next_pos = lposition + direction * 1e-3f;
      }
      pdf += lpdf;
    } else if (light.environment != invalidid) {
      auto  env     = scene.environments[light.environment];
      auto& texture = scene.textures[env.emission_tex];

      // calculate inxex of the texel
      auto wl       = transform_direction(inverse(env.frame), direction);
      auto texcoord = vec2f{
          atan2(wl.z, wl.x) / (2 * pif), acos(clamp(wl.y, -1.0f, 1.0f)) / pif};
      if (texcoord.x < 0) texcoord.x += 1;
      auto i = clamp((int)(texcoord.x * texture.width), 0, texture.width - 1);
      auto j = clamp((int)(texcoord.y * texture.height), 0, texture.height - 1);

      // find the probability of the texel
      /*auto prob = sample_discrete_pdf(
                      light.elements_cdf, j * texture.width + i) /
                  light.elements_cdf.back();*/
      auto prob  = light.probabilities[j * texture.width + i];
      auto angle = (2 * pif / texture.width) * (pif / texture.height) *
                   sin(pif * (j + 0.5f) / texture.height);
      pdf += prob / angle;
    } else {
      pdf += 1 / (4 * pif);
    }
  }

  pdf *= sample_uniform_pdf((int)lights.lights.size());
  return pdf;
}

// Recursive path tracing.
static vec4f shade_pathtrace(const scene_data& scene, const bvh_data& bvh,
    const pathtrace_lights& lights, const ray3f& ray_, rng_state& rng,
    const pathtrace_params& params) {
  // YOUR CODE GOES HERE

  // initialize
  auto radiance = vec3f{0, 0, 0};
  auto weight   = vec3f{1, 1, 1};
  auto ray      = ray_;
  auto hit      = false;

  // trace  path
  for (auto bounce = 0; bounce < max(params.bounces, 4); bounce++) {
    // intersect next point
    auto intersection = intersect_bvh(bvh, scene, ray);
    if (!intersection.hit) {
      radiance += weight * eval_environment(scene, ray.d);
      break;
    }

    // prepare shading point
    auto outgoing = -ray.d;
    auto position = eval_shading_position(scene, intersection, outgoing);
    auto normal   = eval_shading_normal(scene, intersection, outgoing);
    auto material = eval_material(scene, intersection);

    // handle opacity
    if (material.opacity < 1 && rand1f(rng) >= material.opacity) {
      ray = {position + ray.d * 1e-2f, ray.d};
      bounce -= 1;
      continue;
    }

    if (bounce == 0) hit = true;

    // accumulate emission
    radiance += weight * eval_emission(material, normal, outgoing);

    // next direction
    auto incoming = vec3f{0, 0, 0};
    if (!is_delta(material)) {
      incoming = rand1f(rng) > 0.5f
                     ? sample_bsdfcos(
                           material, normal, outgoing, rand1f(rng), rand2f(rng))
                     : sample_lights(scene, lights, position, rand1f(rng),
                           rand1f(rng), rand1f(rng), rand2f(rng));

      if (incoming == vec3f{0, 0, 0}) break;
      weight *= eval_bsdfcos(material, normal, outgoing, incoming) * 2 /
                (sample_bsdfcos_pdf(material, normal, outgoing, incoming) +
                    sample_lights_pdf(scene, bvh, lights, position, incoming));
    } else {
      incoming = sample_delta(material, normal, outgoing, rand1f(rng));
      if (incoming == vec3f{0, 0, 0}) break;
      weight *= eval_delta(material, normal, outgoing, incoming) /
                sample_delta_pdf(material, normal, outgoing, incoming);
    }

    // check weight
    if (weight == vec3f{0, 0, 0} || !isfinite(weight)) break;

    // russian roulette
    if (bounce > 2) {
      auto rr_prob = min((float)0.99, max(weight));
      if (rand1f(rng) >= rr_prob) break;
      weight *= 1 / rr_prob;
    }

    // setup next iteration
    ray = {position, incoming};
  }

  return {radiance.x, radiance.y, radiance.z, hit ? 1.0f : 0.0f};
}

// Recursive path tracing.
static vec4f shade_naive(const scene_data& scene, const bvh_data& bvh,
    const pathtrace_lights& lights, const ray3f& ray_, rng_state& rng,
    const pathtrace_params& params) {
  // YOUR CODE GOES HERE

  // initialize
  auto radiance = vec3f{0, 0, 0};
  auto weight   = vec3f{1, 1, 1};
  auto ray      = ray_;
  auto hit      = false;

  // trace  path
  for (auto bounce = 0; bounce < max(params.bounces, 4); bounce++) {
    // intersect next point
    auto intersection = intersect_bvh(bvh, scene, ray);
    if (!intersection.hit) {
      radiance += weight * eval_environment(scene, ray.d);
      break;
    }

    // prepare shading point
    auto outgoing = -ray.d;
    auto position = eval_shading_position(scene, intersection, outgoing);
    auto normal   = eval_shading_normal(scene, intersection, outgoing);
    auto material = eval_material(scene, intersection);

    // handle opacity
    if (material.opacity < 1 && rand1f(rng) >= material.opacity) {
      ray = {position + ray.d * 1e-2f, ray.d};
      bounce -= 1;
      continue;
    }

    if (bounce == 0) hit = true;

    // accumulate emission
    radiance += weight * eval_emission(material, normal, outgoing);

    // next direction
    auto incoming = vec3f{0, 0, 0};
    if (!is_delta(material)) {
      incoming = sample_bsdfcos(
          material, normal, outgoing, rand1f(rng), rand2f(rng));
      if (incoming == vec3f{0, 0, 0}) break;
      weight *= eval_bsdfcos(material, normal, outgoing, incoming) /
                sample_bsdfcos_pdf(material, normal, outgoing, incoming);
    } else {
      incoming = sample_delta(material, normal, outgoing, rand1f(rng));
      if (incoming == vec3f{0, 0, 0}) break;
      weight *= eval_delta(material, normal, outgoing, incoming) /
                sample_delta_pdf(material, normal, outgoing, incoming);
    }

    // check weight
    if (weight == vec3f{0, 0, 0} || !isfinite(weight)) break;

    // russian roulette
    if (bounce > 2) {
      auto rr_prob = min((float)0.99, max(weight));
      if (rand1f(rng) >= rr_prob) break;
      weight *= 1 / rr_prob;
    }

    // setup next iteration
    ray = {position, incoming};
  }

  return {radiance.x, radiance.y, radiance.z, hit ? 1.0f : 0.0f};
}

// Eyelight for quick previewing.
static vec4f shade_eyelight(const scene_data& scene, const bvh_data& bvh,
    const pathtrace_lights& lights, const ray3f& ray_, rng_state& rng,
    const pathtrace_params& params) {
  // initialize
  auto radiance = vec3f{0, 0, 0};
  auto weight   = vec3f{1, 1, 1};
  auto ray      = ray_;
  auto hit      = false;

  // trace  path
  for (auto bounce = 0; bounce < max(params.bounces, 4); bounce++) {
    // intersect next point
    auto intersection = intersect_bvh(bvh, scene, ray);
    if (!intersection.hit) {
      radiance += weight * eval_environment(scene, ray.d);
      break;
    }

    // prepare shading point
    auto outgoing = -ray.d;
    auto position = eval_shading_position(scene, intersection, outgoing);
    auto normal   = eval_shading_normal(scene, intersection, outgoing);
    auto material = eval_material(scene, intersection);

    // handle opacity
    if (material.opacity < 1 && rand1f(rng) >= material.opacity) {
      ray = {position + ray.d * 1e-2f, ray.d};
      bounce -= 1;
      continue;
    }

    // set hit variables
    if (bounce == 0) hit = true;

    // accumulate emission
    auto incoming = outgoing;
    radiance += weight * eval_emission(material, normal, outgoing);

    // brdf * light
    radiance += weight * pif *
                eval_bsdfcos(material, normal, outgoing, incoming);

    // continue path
    if (!is_delta(material)) break;
    incoming = sample_delta(material, normal, outgoing, rand1f(rng));
    if (incoming == vec3f{0, 0, 0}) break;
    weight *= eval_delta(material, normal, outgoing, incoming) /
              sample_delta_pdf(material, normal, outgoing, incoming);
    if (weight == vec3f{0, 0, 0} || !isfinite(weight)) break;

    // setup next iteration
    ray = {position, incoming};
  }

  return {radiance.x, radiance.y, radiance.z, hit ? 1.0f : 0.0f};
}

// Normal for debugging.
static vec4f shade_normal(const scene_data& scene, const bvh_data& bvh,
    const pathtrace_lights& lights, const ray3f& ray, rng_state& rng,
    const pathtrace_params& params) {
  // intersect next point
  auto intersection = intersect_bvh(bvh, scene, ray);
  if (!intersection.hit) return {0, 0, 0, 0};

  // prepare shading point
  auto outgoing = -ray.d;
  auto normal   = eval_shading_normal(scene, intersection, outgoing);
  return {normal.x, normal.y, normal.z, 1};
}

// Normal for debugging.
static vec4f shade_texcoord(const scene_data& scene, const bvh_data& bvh,
    const pathtrace_lights& lights, const ray3f& ray, rng_state& rng,
    const pathtrace_params& params) {
  // intersect next point
  auto intersection = intersect_bvh(bvh, scene, ray);
  if (!intersection.hit) return {0, 0, 0, 0};

  // prepare shading point
  auto texcoord = eval_texcoord(scene, intersection);
  return {texcoord.x, texcoord.y, 0, 1};
}

// Albedo for denosing
static vec4f shade_albedo(const scene_data& scene, const bvh_data& bvh,
    const pathtrace_lights& lights, const ray3f& ray_, rng_state& rng,
    const pathtrace_params& params) {
  // initialize
  auto radiance = vec3f{0, 0, 0};
  auto ray      = ray_;
  auto hit      = false;

  // trace  path
  for (auto bounce = 0; bounce < min(params.bounces, 4); bounce++) {
    // intersect next point
    auto intersection = intersect_bvh(bvh, scene, ray);
    if (!intersection.hit) {
      radiance += eval_environment(scene, ray.d);
      break;
    }

    // prepare shading point
    auto outgoing = -ray.d;
    auto position = eval_shading_position(scene, intersection, outgoing);
    auto normal   = eval_shading_normal(scene, intersection, outgoing);
    auto material = eval_material(scene, intersection);

    // handle opacity
    if (material.opacity < 1) {
      ray = {position + ray.d * 1e-2f, ray.d};
      bounce -= 1;
      continue;
    }

    if (material.type == material_type::refractive) {
      return {1, 1, 1, 1};
    }

    // set hit variables
    if (bounce == 0) hit = true;

    // accumulate emission
    auto incoming = reflect(outgoing, normal);
    radiance += eval_emission(material, normal, outgoing);

    // brdf * light
    if (material.type == material_type::glossy || material.roughness != 0) radiance += material.color;
    else
      radiance += pif * eval_bsdfcos(material, normal, outgoing, incoming) /
                abs(dot(normal, incoming));

    // continue path
    if (!is_delta(material)) break;

    switch (material.type) {
      case material_type::reflective:
        incoming = sample_reflective(material.color, normal, outgoing);
        break;
      case material_type::transparent: incoming = -outgoing; break;
      default:
        incoming = sample_delta(material, normal, outgoing, rand1f(rng));
        break;
    }

    if (incoming == vec3f{0, 0, 0}) break;

    // setup next iteration
    ray = {position, incoming};
  }

  return {radiance.x, radiance.y, radiance.z, 1};
}

static vec4f shade_color(const scene_data& scene, const bvh_data& bvh,
    const pathtrace_lights& lights, const ray3f& ray, rng_state& rng,
    const pathtrace_params& params) {
  // intersect next point
  auto intersection = intersect_bvh(bvh, scene, ray);
  if (!intersection.hit) return {0, 0, 0, 0};
  // prepare shading point
  auto color = eval_material(scene, intersection).color;
  return {color.x, color.y, color.z, 1};
}

// Trace a single ray from the camera using the given algorithm.
using pathtrace_shader_func = vec4f (*)(const scene_data& scene,
    const bvh_scene& bvh, const pathtrace_lights& lights, const ray3f& ray,
    rng_state& rng, const pathtrace_params& params);
static pathtrace_shader_func get_shader(const pathtrace_params& params) {
  switch (params.shader) {
    case pathtrace_shader_type::pathtrace: return shade_pathtrace;
    case pathtrace_shader_type::naive: return shade_naive;
    case pathtrace_shader_type::eyelight: return shade_eyelight;
    case pathtrace_shader_type::normal: return shade_normal;
    case pathtrace_shader_type::texcoord: return shade_texcoord;
    case pathtrace_shader_type::color: return shade_color;
    case pathtrace_shader_type::albedo: return shade_albedo;
    default: {
      throw std::runtime_error("sampler unknown");
      return nullptr;
    }
  }
}

// Build the bvh acceleration structure.
bvh_scene make_bvh(const scene_data& scene, const pathtrace_params& params) {
  return make_bvh(scene, false, false, params.noparallel);
}

// Init a sequence of random number generators.
pathtrace_state make_state(
    const scene_data& scene, const pathtrace_params& params) {
  auto& camera = scene.cameras[params.camera];
  auto  state  = pathtrace_state{};
  if (camera.aspect >= 1) {
    state.width  = params.resolution;
    state.height = (int)round(params.resolution / camera.aspect);
  } else {
    state.height = params.resolution;
    state.width  = (int)round(params.resolution * camera.aspect);
  }
  state.samples = 0;
  state.image.assign(state.width * state.height, {0, 0, 0, 0});
  state.hits.assign(state.width * state.height, 0);
  state.rngs.assign(state.width * state.height, {});
  auto rng_ = make_rng(1301081);
  for (auto& rng : state.rngs) {
    rng = make_rng(961748941ull, rand1i(rng_, 1 << 31) / 2 + 1);
  }
  return state;
}

// Init trace lights
pathtrace_lights make_lights(
    const scene_data& scene, const pathtrace_params& params) {
  auto lights = pathtrace_lights{};

  // TODO: here we create the cdfs

  for (auto handle = 0; handle < scene.instances.size(); handle++) {
    auto& instance = scene.instances[handle];
    auto& material = scene.materials[instance.material];
    if (material.emission == vec3f{0, 0, 0}) continue;
    auto& shape = scene.shapes[instance.shape];
    if (shape.triangles.empty() && shape.quads.empty()) continue;
    auto& light       = lights.lights.emplace_back();
    light.instance    = handle;
    light.environment = invalidid;
    if (!shape.triangles.empty()) {
      /*light.elements_cdf = vector<float>(shape.triangles.size());*/
      light.probabilities = vector<double>(shape.triangles.size());
      light.alias         = vector<int>(shape.triangles.size());

      for (auto idx = 0; idx < shape.triangles.size(); idx++) {
        auto& t = shape.triangles[idx];
        /*light.elements_cdf[idx] = triangle_area(
            shape.positions[t.x], shape.positions[t.y], shape.positions[t.z]);*/
        auto area = triangle_area(
            shape.positions[t.x], shape.positions[t.y], shape.positions[t.z]);
        light.area += area;
        light.probabilities[idx] = area;

        /*if (idx != 0) light.elements_cdf[idx] += light.elements_cdf[idx -
         * 1];*/
      }
    }
    if (!shape.quads.empty()) {
      // light.elements_cdf = vector<float>(shape.quads.size());

      light.probabilities = vector<double>(shape.quads.size());
      light.alias         = vector<int>(shape.quads.size());

      for (auto idx = 0; idx < shape.quads.size(); idx++) {
        auto& t = shape.quads[idx];
        /*light.elements_cdf[idx] = quad_area(shape.positions[t.x],
            shape.positions[t.y], shape.positions[t.z], shape.positions[t.w]);*/
        auto area = quad_area(shape.positions[t.x], shape.positions[t.y],
            shape.positions[t.z], shape.positions[t.w]);
        light.area += area;
        light.probabilities[idx] = area;

        /*if (idx != 0) light.elements_cdf[idx] += light.elements_cdf[idx -
         * 1];*/
      }
    }

    // normalize probabilities
    for (auto idx = 0; idx < light.probabilities.size(); idx++) {
      light.probabilities[idx] /= light.area;
    }
  }
  for (auto handle = 0; handle < scene.environments.size(); handle++) {
    auto& environment = scene.environments[handle];
    if (environment.emission == vec3f{0, 0, 0}) continue;
    auto& light       = lights.lights.emplace_back();
    light.instance    = invalidid;
    light.environment = handle;
    if (environment.emission_tex != invalidid) {
      auto& texture = scene.textures[environment.emission_tex];
      /*light.elements_cdf = vector<float>(texture.width * texture.height);*/

      auto size = texture.width * texture.height;

      light.probabilities = vector<double>(size);
      light.alias         = vector<int>(size);

      for (auto idx = 0; idx < size; idx++) {
        auto ij    = vec2i{idx % texture.width, idx / texture.width};
        auto th    = (ij.y + 0.5f) * pif / texture.height;
        auto value = lookup_texture(texture, ij.x, ij.y);
        /*light.elements_cdf[idx] = max(value) * sin(th);*/

        light.area += (double)max(value) * sin(th);
        light.probabilities[idx] = (double)max(value) * sin(th);

        /*if (idx != 0) light.elements_cdf[idx] += light.elements_cdf[idx -
         * 1];*/
      }

      // normalize probabilities
      for (auto idx = 0; idx < size; idx++) {
        light.probabilities[idx] /= light.area;
      }
    }
  }

  // init alias method
  for (auto& light : lights.lights) {
    auto average = 1.0 / light.probabilities.size();

    light.probabilities_modified = light.probabilities;
    auto& probs                  = light.probabilities_modified;

    std::deque<int> small;
    std::deque<int> large;

    // populate the stacks with the probabilities
    for (auto i = 0; i < probs.size(); i++) {
      if (probs[i] >= average)
        large.push_back(i);
      else
        small.push_back(i);
    }

    while (!small.empty() && !large.empty()) {
      /* Get the index of the small and the large probabilities. */
      auto less = small.back();
      auto more = large.back();
      small.pop_back();
      large.pop_back();

      /* These probabilities have not yet been scaled up to be such that
       * 1/n is given weight 1.0.  We do this here instead.
       */
      probs[less]       = probs[less] * probs.size();
      light.alias[less] = more;

      /* Decrease the probability of the larger one by the appropriate
       * amount.
       */
      probs[more] = (probs[more] + probs[less]) - average;

      /* If the new probability is less than the average, add it into the
       * small list; otherwise add it to the large list.
       */
      if (probs[more] >= 1.0 / probs.size())
        large.push_back(more);
      else
        small.push_back(more);
    }

    /* At this point, everything is in one list, which means that the
     * remaining probabilities should all be 1/n.  Based on this, set them
     * appropriately.  Due to numerical issues, we can't be sure which
     * stack will hold the entries, so we empty both.
     */
    while (!small.empty()) {
      auto less = small.back();
      small.pop_back();

      probs[less] = 1.0;
    }
    while (!large.empty()) {
      auto more = large.back();
      large.pop_back();
      probs[more] = 1.0;
    }
  }

  // handle progress
  return lights;
}

// Progressively compute an image by calling trace_samples multiple times.
void pathtrace_samples(pathtrace_state& state, const scene_data& scene,
    const bvh_scene& bvh, const pathtrace_lights& lights,
    const pathtrace_params& params) {
  if (state.samples >= params.samples) return;
  auto& camera = scene.cameras[params.camera];
  auto  shader = get_shader(params);
  state.samples += 1;
  if (params.samples == 1) {
    for (auto idx = 0; idx < state.width * state.height; idx++) {
      auto i = idx % state.width, j = idx / state.width;
      auto u = (i + 0.5f) / state.width, v = (j + 0.5f) / state.height;
      auto rng = state.rngs[idx];
      auto ray = eval_camera(
          camera, {u, v}, sample_disk(rand2f(state.rngs[idx])));
      auto radiance = shader(scene, bvh, lights, ray, state.rngs[idx], params);
      if (!isfinite(radiance)) radiance = {0, 0, 0};
      if (params.shader == pathtrace_shader_type::pathtrace ||
          params.shader == pathtrace_shader_type::naive)
        radiance = clamp(radiance, 0.0f, 10.0f);
      state.image[idx] += radiance;
      state.hits[idx] += 1;
    }
  } else if (params.noparallel) {
    for (auto idx = 0; idx < state.width * state.height; idx++) {
      auto i = idx % state.width, j = idx / state.width;
      auto u   = (i + rand1f(state.rngs[idx])) / state.width,
           v   = (j + rand1f(state.rngs[idx])) / state.height;
      auto ray = eval_camera(
          camera, {u, v}, sample_disk(rand2f(state.rngs[idx])));
      auto radiance = shader(scene, bvh, lights, ray, state.rngs[idx], params);
      if (!isfinite(radiance)) radiance = {0, 0, 0};
      if (params.shader == pathtrace_shader_type::pathtrace ||
          params.shader == pathtrace_shader_type::naive)
        radiance = clamp(radiance, 0.0f, 10.0f);
      state.image[idx] += radiance;
      state.hits[idx] += 1;
    }
  } else {
    parallel_for(state.width * state.height, [&](int idx) {
      auto i = idx % state.width, j = idx / state.width;
      auto u   = (i + rand1f(state.rngs[idx])) / state.width,
           v   = (j + rand1f(state.rngs[idx])) / state.height;
      auto ray = eval_camera(
          camera, {u, v}, sample_disk(rand2f(state.rngs[idx])));
      auto radiance = shader(scene, bvh, lights, ray, state.rngs[idx], params);
      if (!isfinite(radiance)) radiance = {0, 0, 0};
      if (params.shader == pathtrace_shader_type::pathtrace ||
          params.shader == pathtrace_shader_type::naive)
        radiance = clamp(radiance, 0.0f, 10.0f);
      state.image[idx] += radiance;
      state.hits[idx] += 1;
    });
  }
}

// Check image type
static void check_image(
    const color_image& image, int width, int height, bool linear) {
  if (image.width != width || image.height != height)
    throw std::invalid_argument{"image should have the same size"};
  if (image.linear != linear)
    throw std::invalid_argument{
        linear ? "expected linear image" : "expected srgb image"};
}

// Get resulting render
color_image get_render(const pathtrace_state& state) {
  auto image = make_image(state.width, state.height, true);
  get_render(image, state);
  return image;
}
void get_render(color_image& image, const pathtrace_state& state) {
  check_image(image, state.width, state.height, true);
  auto scale = 1.0f / (float)state.samples;
  for (auto idx = 0; idx < state.width * state.height; idx++) {
    image.pixels[idx] = state.image[idx] * scale;
  }
}

}  // namespace yocto
