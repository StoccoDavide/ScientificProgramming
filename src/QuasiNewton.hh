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

  typedef Eigen::DiagonalMatrix<real, 3>           scale;     //!< 3D scaling transformation type
  typedef Eigen::Translation<real, 3>              translate; //!< 3D translation transformation type
  typedef Eigen::AngleAxis<real>                   angleaxis; //!< 3D rotation transformation type
  typedef Eigen::Transform<real, 3, Eigen::Affine> affine;    //!< 3D affine transformation type

} // namespace QuasiNewton

#include "QuasiNewton/Solver.hxx"

#endif

///
/// eof: QuasiNewton.hh
///
