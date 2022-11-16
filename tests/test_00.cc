#include "QuasiNewton.hh"

using namespace QuasiNewton;

int
main(void)
{
  // Eirola-Nevanlinna Method
  EirolaNevanlinna<vec2, mat2> test_ENM_2;
  EirolaNevanlinna<vec3, mat3> test_ENM_3;
  EirolaNevanlinna<vec4, mat4> test_ENM_4;

  // Broyden's Ugly Method
  BroydenUgly<vec2, mat2> test_BUM_2;
  BroydenUgly<vec3, mat3> test_BUM_3;
  BroydenUgly<vec4, mat4> test_BUM_4;

  // Broyden's Bad Method
  BroydenBad<vec2, mat2> test_BBM_2;
  BroydenBad<vec3, mat3> test_BBM_3;
  BroydenBad<vec4, mat4> test_BBM_4;

  // Broyden's Good Method
  BroydenGood<vec2, mat2> test_BGM_2;
  BroydenGood<vec3, mat3> test_BGM_3;
  BroydenGood<vec4, mat4> test_BGM_4;

  // Broyden's Good Method
  BroydenCombined<vec2, mat2> test_BCM_2;
  BroydenCombined<vec3, mat3> test_BCM_3;
  BroydenCombined<vec4, mat4> test_BCM_4;

  // Greenstadt's 1st Method
  Greenstadt1<vec2, mat2> test_G1M_2;
  Greenstadt1<vec3, mat3> test_G1M_3;
  Greenstadt1<vec4, mat4> test_G1M_4;

  // Greenstadt's 2nd Method
  Greenstadt2<vec2, mat2> test_G2M_2;
  Greenstadt2<vec3, mat3> test_G2M_3;
  Greenstadt2<vec4, mat4> test_G2M_4;

  return 0;
}
