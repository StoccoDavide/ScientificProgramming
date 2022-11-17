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
/// file: QuasiNewton.hh
///

#pragma once

#ifndef INCLUDE_QUASINEWTON
#define INCLUDE_QUASINEWTON

// Print QuasiNewton errors
#ifndef QUASINEWTON_ERROR
#define QUASINEWTON_ERROR(MSG)                \
  {                                     \
    std::ostringstream os;              \
    os << MSG;                          \
    throw std::runtime_error(os.str()); \
  }
#endif

// Assert for QuasiNewton
#ifndef QUASINEWTON_ASSERT
#define QUASINEWTON_ASSERT(COND, MSG) \
  if (!(COND))                  \
  QUASINEWTON_ERROR(MSG)
#endif

// Standard libraries
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>

// Eigen library
#include <Eigen/Dense>

//! Namespace containing all Quasi-Newton typedefs, classes and routines
namespace QuasiNewton
{

  /*\
   |   _____                     _       __
   |  |_   _|   _ _ __   ___  __| | ___ / _|___
   |    | || | | | '_ \ / _ \/ _` |/ _ \ |_/ __|
   |    | || |_| | |_) |  __/ (_| |  __/  _\__ \
   |    |_| \__, | .__/ \___|\__,_|\___|_| |___/
   |        |___/|_|
  \*/

  typedef double       real;       //!< Real number type
  typedef int          integer;    //!< Integer number type
  typedef std::ostream out_stream; //!< Output stream type

  typedef Eigen::Matrix<real, 2, 1>                           vec2; //!< 2x1 vector of real number type (column vector)
  typedef Eigen::Matrix<real, 2, 2>                           mat2; //!< 2x2 matrix of real number type
  typedef Eigen::Matrix<real, 3, 1>                           vec3; //!< 3x1 vector of real number type (column vector)
  typedef Eigen::Matrix<real, 3, 3>                           mat3; //!< 3x3 matrix of real number type
  typedef Eigen::Matrix<real, 4, 1>                           vec4; //!< 4x1 vector of real number type (column vector)
  typedef Eigen::Matrix<real, 4, 4>                           mat4; //!< 4x4 matrix of real number type
  typedef Eigen::Matrix<real, 5, 1>                           vec5; //!< 5x1 vector of real number type (column vector)
  typedef Eigen::Matrix<real, 5, 5>                           mat5; //!< 5x5 matrix of real number type
  typedef Eigen::Matrix<real, 6, 1>                           vec6; //!< 6x1 vector of real number type (column vector)
  typedef Eigen::Matrix<real, 6, 6>                           mat6; //!< 6x6 matrix of real number type
  typedef Eigen::Matrix<real, 7, 1>                           vec7; //!< 7x1 vector of real number type (column vector)
  typedef Eigen::Matrix<real, 7, 7>                           mat7; //!< 7x7 matrix of real number type
  typedef Eigen::Matrix<real, 8, 1>                           vec8; //!< 8x1 vector of real number type (column vector)
  typedef Eigen::Matrix<real, 8, 8>                           mat8; //!< 8x8 matrix of real number type
  typedef Eigen::Matrix<real, 9, 1>                           vec9; //!< 9x1 vector of real number type (column vector)
  typedef Eigen::Matrix<real, 9, 9>                           mat9; //!< 9x9 matrix of real number type
  typedef Eigen::Matrix<real, 10, 1>                          vec10; //!< 10x1 vector of real number type (column vector)
  typedef Eigen::Matrix<real, 10, 10>                         mat10; //!< 10x10 matrix of real number type
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1>              vecN; //!< Nx1 vector of real number type (column vector)
  typedef Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> matN; //!< NxN matrix of real number type

  /*\
   |    ____                _              _
   |   / ___|___  _ __  ___| |_ __ _ _ __ | |_ ___
   |  | |   / _ \| '_ \/ __| __/ _` | '_ \| __/ __|
   |  | |__| (_) | | | \__ \ || (_| | | | | |_\__ \
   |   \____\___/|_| |_|___/\__\__,_|_| |_|\__|___/
   |
  \*/

  static real const EPSILON_MACHINE = std::numeric_limits<real>::epsilon();      //!< Machine epsilon epsilon static constant value
  static real const EPSILON_HIGH    = 1.0E-16;                                   //!< High precision epsilon static constant value
  static real const EPSILON_MEDIUM  = 1.0E-10;                                   //!< Medium precision epsilon static constant value
  static real const EPSILON_LOW     = 1.0E-07;                                   //!< Low precision epsilon static constant value
  static real const EPSILON         = EPSILON_MEDIUM;                            //!< Standard precision epsilon static constant value
  static real const INFTY           = std::numeric_limits<real>::infinity();     //!< Infinity static constant value
  static real const QUIET_NAN       = std::numeric_limits<real>::quiet_NaN();    //!< Not-a-Number static constant value
  static real const PI              = real(3.141592653589793238462643383279500); //!< Pi static constant value
  static real const PIDIV180        = real(0.017453292519943295769236907684886); //!< Pi/180 static constant value

  static vec2 const UNITX_VEC2    = vec2::UnitX();             //!< X axis unit vec2 static constant object
  static vec2 const UNITY_VEC2    = vec2::UnitY();             //!< Y axis unit vec2 static constant object
  static vec2 const NAN_VEC2      = vec2::Constant(QUIET_NAN); //!< Not-a-Number vec2 static constant object
  static mat2 const NAN_MAT2      = mat2::Constant(QUIET_NAN); //!< Not-a-Number mat2 static constant object
  static vec2 const ZEROS_VEC2    = vec2::Constant(0.0);       //!< Zeros vec2 static constant object
  static mat2 const ZEROS_MAT2    = mat2::Constant(0.0);       //!< Zeros mat2 static constant object
  static vec2 const ONES_VEC2     = vec2::Constant(1.0);       //!< Ones vec2 static constant object
  static mat2 const ONES_MAT2     = mat2::Constant(1.0);       //!< Ones mat2 static constant object
  static mat2 const IDENTITY_MAT2 = mat2::Identity();          //!< Identity mat2 static constant object

  static vec3 const UNITX_VEC3    = vec3::UnitX();             //!< X axis unit vec3 type
  static vec3 const UNITY_VEC3    = vec3::UnitY();             //!< Y axis unit vec3 type
  static vec3 const UNITZ_VEC3    = vec3::UnitZ();             //!< Z axis unit vec3 type
  static vec3 const NAN_VEC3      = vec3::Constant(QUIET_NAN); //!< Not-a-Number vec3 type
  static mat3 const NAN_MAT3      = mat3::Constant(QUIET_NAN); //!< Not-a-Number mat3 type
  static vec3 const ZEROS_VEC3    = vec3::Constant(0.0);       //!< Zeros vec3 type
  static mat3 const ZEROS_MAT3    = mat3::Constant(0.0);       //!< Zeros mat3 type
  static vec3 const ONES_VEC3     = vec3::Constant(1.0);       //!< Ones vec3 type
  static mat3 const ONES_MAT3     = mat3::Constant(1.0);       //!< Ones mat3 type
  static mat3 const IDENTITY_MAT3 = mat3::Identity();          //!< Identity mat3 type

  static vec4 const UNITX_VEC4    = vec4::UnitX();             //!< X axis unit vec4 type
  static vec4 const UNITY_VEC4    = vec4::UnitY();             //!< Y axis unit vec4 type
  static vec4 const UNITZ_VEC4    = vec4::UnitZ();             //!< Z axis unit vec4 type
  static vec4 const UNITW_VEC4    = vec4::UnitW();             //!< W axis unit vec4 type
  static vec4 const NAN_VEC4      = vec4::Constant(QUIET_NAN); //!< Not-a-Number vec4 type
  static mat4 const NAN_MAT4      = mat4::Constant(QUIET_NAN); //!< Not-a-Number mat4 type
  static vec4 const ZEROS_VEC4    = vec4::Constant(0.0);       //!< Zeros vec4 type
  static mat4 const ZEROS_MAT4    = mat4::Constant(0.0);       //!< Zeros mat4 type
  static vec4 const ONES_VEC4     = vec4::Constant(1.0);       //!< Ones vec4 type
  static mat4 const ONES_MAT4     = mat4::Constant(1.0);       //!< Ones mat4 type
  static mat4 const IDENTITY_MAT4 = mat4::Identity();          //!< Identity mat4 type

  static real DUMMY_REAL = QUIET_NAN; //!< Dummy vec2 type static non-const object
  static vec2 DUMMY_VEC2 = NAN_VEC2;  //!< Dummy vec2 type static non-const object
  static vec3 DUMMY_VEC3 = NAN_VEC3;  //!< Dummy vec3 type static non-const object
  static vec4 DUMMY_VEC4 = NAN_VEC4;  //!< Dummy vec4 type static non-const object
  static mat2 DUMMY_MAT2 = NAN_MAT2;  //!< Dummy mat2 type static non-const object
  static mat3 DUMMY_MAT3 = NAN_MAT3;  //!< Dummy mat3 type static non-const object
  static mat4 DUMMY_MAT4 = NAN_MAT4;  //!< Dummy mat4 type static non-const object

} // namespace QuasiNewton

#include "QuasiNewton/Solver.hxx"

#endif

///
/// eof: QuasiNewton.hh
///
