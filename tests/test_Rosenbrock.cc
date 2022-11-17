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
  vec2 Xini(-1.75, 3.0);
  mat2 Gini(IDENTITY_MAT2);
  vec2 Xout(QUIET_NAN, QUIET_NAN);

  // Solve without damping
  std::cout
      << "SOLVER WITHOUT DAMPING" << std::endl
      << std::endl;
  for (size_t i = 0; i < SolversPtrVec.size(); ++i)
  {
    SolversPtrVec[i]->solve(
      [](vec2 const & X, vec2 & F)
        {
          real a = real(1.0);
          real b = real(100.0);
          F(0) = a*(1 - X(0));
          F(1) = b*(X(1) - X(0)*X(0));
        },
      Xini,
      Gini,
      Xout
    );
    std::cout
      << "SOLVER "<< i << ": " << SolversPtrVec[i]->method() << std::endl
      << "\tConverged   = " << SolversPtrVec[i]->outConverged() << " " << std::endl
      << "\tSolution    = [ " << Xout.transpose() << " ]'" << std::endl
      << "\tEvaluations = " << SolversPtrVec[i]->outEvaluations() << std::endl
      << "\tIterations  = " << SolversPtrVec[i]->outIterations() << std::endl
      << "\tRelaxations = " << SolversPtrVec[i]->outRelaxations() << std::endl
      << std::endl;
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
    std::cout
      << "SOLVER "<< i << ": " << SolversPtrVec[i]->method() << std::endl
      << "\tConverged   = " << SolversPtrVec[i]->outConverged() << " " << std::endl
      << "\tSolution    = [ " << Xout.transpose() << " ]'" << std::endl
      << "\tEvaluations = " << SolversPtrVec[i]->outEvaluations() << std::endl
      << "\tIterations  = " << SolversPtrVec[i]->outIterations() << std::endl
      << "\tRelaxations = " << SolversPtrVec[i]->outRelaxations() << std::endl
      << std::endl;
  }

  return 0;
}
