/*
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * MIT License                                                                     *
 *                                                                                 *
 * Copyright (c) 2022, Davide Stocco                                               *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy    *
 * of this software and associated documentation files (the "Software"), to deal   *
 * in the Software without restriction, including without limitation the rights    *
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell       *
 * copies of the Software, and to permit persons to whom the Software is           *
 * furnished to do so, subject to the following conditions:                        *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included in all  *
 * copies or substantial portions of the Software.                                 *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR      *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,   *
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE   *
 * SOFTWARE.                                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

///
/// file: QuasiNewton.cc
///

#include "QuasiNewton.hh"

namespace QuasiNewton
{
  // Eirola-Nevanlinna Method
  template class EirolaNevanlinna<vec2, mat2>;
  template class EirolaNevanlinna<vec3, mat3>;
  template class EirolaNevanlinna<vec4, mat4>;
  template class EirolaNevanlinna<vec5, mat5>;
  template class EirolaNevanlinna<vec6, mat6>;
  template class EirolaNevanlinna<vec7, mat7>;
  template class EirolaNevanlinna<vec8, mat8>;
  template class EirolaNevanlinna<vec9, mat9>;
  template class EirolaNevanlinna<vec10, mat10>;

  // Broyden's Ugly Method
  template class BroydenUgly<vec2, mat2>;
  template class BroydenUgly<vec3, mat3>;
  template class BroydenUgly<vec4, mat4>;
  template class BroydenUgly<vec5, mat5>;
  template class BroydenUgly<vec6, mat6>;
  template class BroydenUgly<vec7, mat7>;
  template class BroydenUgly<vec8, mat8>;
  template class BroydenUgly<vec9, mat9>;
  template class BroydenUgly<vec10, mat10>;

  // Broyden's Bad Method
  template class BroydenBad<vec2, mat2>;
  template class BroydenBad<vec3, mat3>;
  template class BroydenBad<vec4, mat4>;
  template class BroydenBad<vec5, mat5>;
  template class BroydenBad<vec6, mat6>;
  template class BroydenBad<vec7, mat7>;
  template class BroydenBad<vec8, mat8>;
  template class BroydenBad<vec9, mat9>;
  template class BroydenBad<vec10, mat10>;

  // Broyden's Good Method
  template class BroydenGood<vec2, mat2>;
  template class BroydenGood<vec3, mat3>;
  template class BroydenGood<vec4, mat4>;
  template class BroydenGood<vec5, mat5>;
  template class BroydenGood<vec6, mat6>;
  template class BroydenGood<vec7, mat7>;
  template class BroydenGood<vec8, mat8>;
  template class BroydenGood<vec9, mat9>;
  template class BroydenGood<vec10, mat10>;

  // Broyden's Good Method
  template class BroydenCombined<vec2, mat2>;
  template class BroydenCombined<vec3, mat3>;
  template class BroydenCombined<vec4, mat4>;
  template class BroydenCombined<vec5, mat5>;
  template class BroydenCombined<vec6, mat6>;
  template class BroydenCombined<vec7, mat7>;
  template class BroydenCombined<vec8, mat8>;
  template class BroydenCombined<vec9, mat9>;
  template class BroydenCombined<vec10, mat10>;

  // Greenstadt's 1st Method
  template class Greenstadt1<vec2, mat2>;
  template class Greenstadt1<vec3, mat3>;
  template class Greenstadt1<vec4, mat4>;
  template class Greenstadt1<vec5, mat5>;
  template class Greenstadt1<vec6, mat6>;
  template class Greenstadt1<vec7, mat7>;
  template class Greenstadt1<vec8, mat8>;
  template class Greenstadt1<vec9, mat9>;
  template class Greenstadt1<vec10, mat10>;

  // Greenstadt's 2nd Method
  template class Greenstadt2<vec2, mat2>;
  template class Greenstadt2<vec3, mat3>;
  template class Greenstadt2<vec4, mat4>;
  template class Greenstadt2<vec5, mat5>;
  template class Greenstadt2<vec6, mat6>;
  template class Greenstadt2<vec7, mat7>;
  template class Greenstadt2<vec8, mat8>;
  template class Greenstadt2<vec9, mat9>;
  template class Greenstadt2<vec10, mat10>;


} // namespace QuasiNewton


///
/// eof: QuasiNewton.cc
///
