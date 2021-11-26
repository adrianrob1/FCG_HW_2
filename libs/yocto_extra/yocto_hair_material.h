#ifndef _YOCTO_HAIR_MATERIAL_H_
#define _YOCTO_HAIR_MATERIAL_H_

// -----------------------------------------------------------------------------
// INCLUDE
// -----------------------------------------------------------------------------

#include <yocto/yocto_math.h>

// -----------------------------------------------------------------------------
// MATERIALS
// -----------------------------------------------------------------------------

namespace yocto::hair {

inline const int p_max = 3;

struct hair_data {
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

#endif