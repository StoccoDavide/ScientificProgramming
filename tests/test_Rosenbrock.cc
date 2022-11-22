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
/// file: test_Rosenbrock.cc
///

#include "QuasiNewton.hh"

using namespace QuasiNewton;

int
main(void)
{
  // Solvers pointer vector
  std::vector<Solver<vec2, mat2> *> SolversPtrVec;
  SolversPtrVec.push_back( new EirolaNevanlinna<vec2, mat2>() );
  SolversPtrVec.push_back( new BroydenUgly<vec2, mat2>()      );
  SolversPtrVec.push_back( new BroydenBad<vec2, mat2>()       );
  SolversPtrVec.push_back( new BroydenGood<vec2, mat2>()      );
  SolversPtrVec.push_back( new BroydenCombined<vec2, mat2>()  );
  SolversPtrVec.push_back( new Greenstadt1<vec2, mat2>()      );
  SolversPtrVec.push_back( new Greenstadt2<vec2, mat2>()      );

  // Starting entries
  vec2 Xini(1.5, 2.0);
  mat2 Gini(IDENTITY_MAT2);
  vec2 Xout(QUIET_NAN, QUIET_NAN);
  real a = real(1.0);
  real b = real(10.0);
  //Gini << -a,                   real(0.0),
  //        -real(2.0)*b*Xini(0), b;
  Gini << -a, real(0.0),
          real(0.0), b;
  Gini = (Gini.inverse()).eval();

  // Solve without damping
  std::cout
      << "SOLVER WITHOUT DAMPING" << std::endl
      << std::endl;
  for (size_t i = 0; i < SolversPtrVec.size(); ++i)
  {
    SolversPtrVec[i]->solve(
      [a, b](vec2 const & X, vec2 & F)
        {
          F(0) = a*(1 - X(0));
          F(1) = b*(X(1) - X(0)*X(0));
        },
      Xini,
      Gini,
      Xout
    );

    // Display the results on command line
    std::cout
      << "SOLVER "<< i << ": " << SolversPtrVec[i]->method() << std::endl
      << "\tConverged   = " << SolversPtrVec[i]->outConverged() << " " << std::endl
      << "\tSolution    = [ " << Xout.transpose() << " ]'" << std::endl
      << "\tEvaluations = " << SolversPtrVec[i]->outEvaluations() << std::endl
      << "\tIterations  = " << SolversPtrVec[i]->outIterations() << std::endl
      << "\tRelaxations = " << SolversPtrVec[i]->outRelaxations() << std::endl
      << std::endl;

    // Print the results on file
    std::ofstream file;
    file.open("./tests/Rosenbrock/Rosenbrock_" + SolversPtrVec[i]->method() + ".csv");
    file
      << std::scientific
      << std::setprecision(6)
      << a << ", " << b << std::endl;
    for (integer k = 0; k < SolversPtrVec[i]->outIterations(); ++k)
    {
      for (integer j = 0; j < 2; ++j)
      {
      file
        << SolversPtrVec[i]->outHistoryX()[k](j) << ", "
        << SolversPtrVec[i]->outHistoryF()[k](j) << ", "
        << SolversPtrVec[i]->outHistoryD()[k](j) << std::endl;
      }
    }
    file.close();
  }

  // Solve with damping
  std::cout
      << "SOLVER WITH DAMPING" << std::endl
      << std::endl;
  for (size_t i = 0; i < SolversPtrVec.size(); ++i)
  {
    SolversPtrVec[i]->solveDumped(
      [](vec2 const & X, vec2 & F)
        {
          real a = real(1.0);
          real b = real(10.0);
          F(0) = a*(1 - X(0));
          F(1) = b*(X(1) - X(0)*X(0));
        },
      Xini,
      Gini,
      Xout
    );

    // Display the results on command line
    std::cout
      << "SOLVER "<< i << ": " << SolversPtrVec[i]->method() << std::endl
      << "\tConverged   = " << SolversPtrVec[i]->outConverged() << " " << std::endl
      << "\tSolution    = [ " << Xout.transpose() << " ]'" << std::endl
      << "\tEvaluations = " << SolversPtrVec[i]->outEvaluations() << std::endl
      << "\tIterations  = " << SolversPtrVec[i]->outIterations() << std::endl
      << "\tRelaxations = " << SolversPtrVec[i]->outRelaxations() << std::endl
      << std::endl;

    // Print the results on file
    std::ofstream file;
    file.open("./tests/Rosenbrock/Rosenbrock_Dumped_" + SolversPtrVec[i]->method() + ".csv");
    file
      << std::scientific
      << std::setprecision(6)
      << a << ", " << b << std::endl;
    for (integer k = 0; k < SolversPtrVec[i]->outIterations(); ++k)
    {
      for (integer j = 0; j < 2; ++j)
      {
      file
        << SolversPtrVec[i]->outHistoryX()[k](j) << ", "
        << SolversPtrVec[i]->outHistoryF()[k](j) << ", "
        << SolversPtrVec[i]->outHistoryD()[k](j) << std::endl;
      }
    }
    file.close();
  }

  return 0;
}

///
/// eof: test_Rosenbrock.cc
///
