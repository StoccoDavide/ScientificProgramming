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
/// file: Solver.hxx
///

#pragma once

#ifndef INCLUDE_QUASINEWTON_SOLVER
#define INCLUDE_QUASINEWTON_SOLVER

namespace QuasiNewton
{

  /*\
   |   ____        _                ____                 _____                 _   _
   |  / ___|  ___ | |_   _____ _ __| __ )  __ _ ___  ___|  ___|   _ _ __   ___| |_(_) ___  _ __
   |  \___ \ / _ \| \ \ / / _ \ '__|  _ \ / _` / __|/ _ \ |_ | | | | '_ \ / __| __| |/ _ \| '_ \
   |   ___) | (_) | |\ V /  __/ |  | |_) | (_| \__ \  __/  _|| |_| | | | | (__| |_| | (_) | | | |
   |  |____/ \___/|_| \_/ \___|_|  |____/ \__,_|___/\___|_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|
   |
  \*/

  //! QuasiNewton solver base function class container
  template <typename vecT>
  class SolverBaseFunction
  {
  public:
    //! Evaluate function F(X)
    virtual
    void
    eval(
      vecT const & X, //!< Input point
      vecT       & F  //!< Output value
    ) const = 0;

    //! Evaluate function F(X) through operator()
    void
    operator()(
      vecT const & X, //!< Input point
      vecT       & F  //!< Output value
    )
    const
    {
      return this->eval(X,F);
    }

  }; // class SolverBaseFunction

  /*\
   |   ____        _                _____                 _   _
   |  / ___|  ___ | |_   _____ _ __|  ___|   _ _ __   ___| |_(_) ___  _ __
   |  \___ \ / _ \| \ \ / / _ \ '__| |_ | | | | '_ \ / __| __| |/ _ \| '_ \
   |   ___) | (_) | |\ V /  __/ |  |  _|| |_| | | | | (__| |_| | (_) | | | |
   |  |____/ \___/|_| \_/ \___|_|  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|
   |
  \*/

  //! QuasiNewton solver function class container
  template <typename vecT, typename ptrT>
  class SolverFunction : public SolverBaseFunction<vecT>
  {
    ptrT m_FunctionPtr;
  public:
    //! Class contructor
    SolverFunction(
      ptrT t_FunctionPtr
    )
      : m_FunctionPtr(t_FunctionPtr)
    {
    }

    //! Evaluate function F(X)
    void
    eval(
      vecT const & X, //!< Input point
      vecT       & F  //!< Output value
    )
    const override
    {
      this->m_FunctionPtr(X, F);
    }

  }; // class SolverFunction

  /*\
   |   ____        _
   |  / ___|  ___ | |_   _____ _ __
   |  \___ \ / _ \| \ \ / / _ \ '__|
   |   ___) | (_) | |\ V /  __/ |
   |  |____/ \___/|_| \_/ \___|_|
   |
  \*/

  //! QuasiNewton solver class container
  template <typename vecT, typename matT>
  class Solver
  {
  private:
    SolverBaseFunction<vecT> *
            m_FunctionPtr    = nullptr;      //!< Function pointer
    real    m_Tolerance      = real(1.0E-8); //!< Algorithm tolerance
    integer m_MaxIterations  = integer(100); //!< Maximum number of allowed algorithm iterations
    integer m_MaxEvaluations = integer(100); //!< Maximum number of allowed function evaluations
    integer m_MaxRelaxations = integer(10);  //!< Maximum number of allowed algorithm relaxations
    bool    m_Verbose        = false;        //!< Verbose mode boolean
    integer m_Iterations     = integer(0);   //!< Algorithm iterations number
    integer m_Evaluations    = integer(0);   //!< Function evaluations number
    integer m_Relaxations    = integer(0);   //!< Algorithm relaxations number
    real    m_Residuals      = real(0.0);    //!< Function residuals
    bool    m_Converged      = false;        //!< Convergence boolean
    real    m_Alpha          = real(0.9);    //!< Relaxation factor

    std::vector<vecT> m_HistoryX; //!< History for X values
    std::vector<vecT> m_HistoryF; //!< History for F values
    std::vector<vecT> m_HistoryD; //!< History for D values

  public:
    //! Class constructor
    Solver(void)
    {
      integer size = this->m_MaxIterations * this->m_MaxRelaxations;
      this->m_HistoryX.reserve(size);
      this->m_HistoryF.reserve(size);
      this->m_HistoryD.reserve(size);
    }

    //! Set algorithm tolerance
    void
    setTolerance(
      real t_Tolerance //!< The algorithm tolerance
    )
    {
      #define CMD "QuasiNewton::Solver::setTolerance(...): "

      QUASINEWTON_ASSERT(
        !std::isnan(t_Tolerance) && std::isfinite(t_Tolerance) && t_Tolerance > real(0.0),
        CMD "invalid input detected.");

      this->m_Tolerance = t_Tolerance;

      #undef CMD
    }

    //!< Get algorithm tolerance
    real
    getTolerance(void)
    const
    {
      return this->m_Tolerance;
    }

    //! Set maximum number of allowed algorithm iterations
    void
    setMaxIterations(
      integer t_MaxIterations //! The maximum number of allowed algorithm iterations
    )
    {
      #define CMD "QuasiNewton::Solver::setMaxIterations(...): "

      QUASINEWTON_ASSERT(
        !std::isnan(t_MaxIterations) && std::isfinite(t_MaxIterations) && t_MaxIterations > integer(0),
        CMD "invalid input detected.");

      this->m_MaxIterations = t_MaxIterations;

      #undef CMD
    }

    //! Get maximum number of allowed algorithm iterations
    integer
    getMaxIterations(void)
    const
    {
      return this->m_MaxEvaluations;
    }

    //! Set maximum number of allowed function evaluations
    void
    setMaxEvaluations(
      integer t_MaxEvaluations //!< The maximum number of allowed function evaluations
    )
    {
      #define CMD "QuasiNewton::Solver::setMaxEvaluations(...): "

      QUASINEWTON_ASSERT(
        !std::isnan(t_MaxEvaluations) && std::isfinite(t_MaxEvaluations) && t_MaxEvaluations > integer(0),
        CMD "invalid input detected.");

      this->m_MaxEvaluations = t_MaxEvaluations;

      #undef CMD
    }

    //! Get maximum number of allowed function evaluations
    integer
    getMaxEvaluations(void)
    const
    {
      return this->m_MaxEvaluations;
    }

    //! Set maximum number of allowed algorithm relaxations
    void
    setMaxRelaxations(
      integer t_MaxRelaxations //!< The maximum number of allowed algorithm relaxations
    )
    {
      #define CMD "QuasiNewton::Solver::setMaxRelaxations(...): "

      QUASINEWTON_ASSERT(
        !std::isnan(t_MaxRelaxations) && std::isfinite(t_MaxRelaxations) && t_MaxRelaxations > integer(0),
        CMD "invalid input detected.");

      this->m_MaxRelaxations = t_MaxRelaxations;

      #undef CMD
    }

    //! Get maximum number of allowed algorithm relaxations
    integer
    getMaxRelaxations(void)
    const
    {
      return this->m_MaxRelaxations;
    }

    //! Set relaxation factor
    void
    setAlpha(
      real t_Alpha //!< The relaxation factor
    )
    {
      #define CMD "QuasiNewton::Solver::setAlpha(...): "

      QUASINEWTON_ASSERT(
        !std::isnan(t_Alpha) && std::isfinite(t_Alpha) &&  real(0.0) <= t_Alpha && t_Alpha <= real(1.0),
        CMD "invalid input detected.");

      this->m_Alpha = t_Alpha;

      #undef CMD
    }

    //! Get relaxation factor
    real
    getAlpha(void)
    const
    {
      return this->m_Alpha;
    }

    //! Enable verbose mode
    void
    enableVerbose(void)
    {
      this->m_Verbose = true;
    }

    //! Disable verbose mode
    void
    disableVerbose(void)
    {
      this->m_Verbose = false;
    }

    //! Get algorithm iterations number
    integer
    outIterations(void)
    const
    {
      return this->m_Iterations;
    }

    //! Get function evaluations number
    integer
    outEvaluations(void)
    const
    {
      return this->m_Evaluations;
    }

    //! Get algorithm relaxations number
    integer
    outRelaxations(void)
    const
    {
      return this->m_Relaxations;
    }

    //! Get function evaluations number
    real
    outResiduals(void)
    const
    {
      return this->m_Residuals;
    }

    //! Get convergence boolean value
    bool
    outConverged(void)
    const
    {
      return this->m_Converged;
    }

    //! Get history of X values
    std::vector<vecT> const &
    outHistoryX(void)
    const
    {
      return this->m_HistoryX;
    }

    //! Get history of F values
    std::vector<vecT> const &
    outHistoryF(void)
    const
    {
      return this->m_HistoryF;
    }

    //! Get history of D values
    std::vector<vecT> const &
    outHistoryD(void)
    const
    {
      return this->m_HistoryD;
    }

    //! Solve non-linear system of equations F(X)=0
    bool
    solve(
      SolverBaseFunction<vecT> * t_FunctionPtr, //!< The fuction pointer
      vecT               const & Xini,          //!< Initialization point
      matT               const & Jini,          //!< Initialization jacobian
      vecT                     & Xout           //!< Output point
    )
    {
      this->m_FunctionPtr = t_FunctionPtr;
      return this->solve(Xini, Jini, Xout);
    }

    //! Solve non-linear system of equations F(X)=0 with dumping
    bool
    solveDumped(
      SolverBaseFunction<vecT> * t_FunctionPtr, //!< The fuction pointer
      vecT               const & Xini,          //!< Initialization point
      matT               const & Jini,          //!< Initialization jacobian
      vecT                     & Xout           //!< Output point
    )
    {
      this->m_FunctionPtr = t_FunctionPtr;
      return this->solveDumped(Xini, Jini, Xout);
    }

    //! Solve non-linear system of equations F(X)=0
    template <typename ptrT>
    bool
    solve(
      ptrT         FunctionPtr, //!< The fuction pointer
      vecT const & Xini,        //!< Initialization point
      matT const & Jini,        //!< Initialization jacobian approximation
      vecT       & Xout         //!< Output point
    )
    {
      SolverFunction<vecT,ptrT> tmp(FunctionPtr);
      this->m_FunctionPtr = &tmp;
      return this->solve(Xini, Jini, Xout);
    }

    //! Solve non-linear system of equations F(X)=0 with dumping
    template <typename ptrT>
    bool
    solveDumped(
      ptrT         FunctionPtr, //!< The fuction pointer
      vecT const & Xini,        //!< Initialization point
      matT const & Jini,        //!< Initialization jacobian approximation
      vecT       & Xout         //!< Output point
    )
    {
      SolverFunction<vecT,ptrT> tmp(FunctionPtr);
      this->m_FunctionPtr = &tmp;
      return this->solveDumped(Xini, Jini, Xout);
    }

    //! Return method name
    virtual
    std::string
    method(void) const = 0;

  protected:
    //! Reset solver internal counter and variables
    void
    reset(void)
    {
      this->m_Iterations  = integer(0);
      this->m_Evaluations = integer(0);
      this->m_Relaxations = integer(0);
      this->m_Residuals   = real(0.0);
      this->m_Converged   = false;
      this->m_HistoryX.clear();
      this->m_HistoryF.clear();
      this->m_HistoryD.clear();
    }

    //! Perform function evaluation
    void
    eval(
      vecT const & X, //!< Input point
      vecT       & F  //!< Output function
    )
    {
      ++this->m_Evaluations;
      this->m_FunctionPtr->eval(X,F);
    }

    //! Solve non-linear system of equations F(X)=0
    bool
    solve(
      vecT const & Xini, //!< Initialization point
      matT const & Gini, //!< Initialization jacobian approximation
      vecT       & Xout  //!< Output point
    )
    {
      // Setup internal variables
      this->reset();

      // Initialize variables
      matT G0, G1;
      vecT X0, F0, D0, DX0, DF0, X1, F1, D1, DX1, DF1;
      real F0_norm = real(0.0);
      real D0_norm = real(0.0);

      // Set initial iteration
      X0 = Xini;
      this->eval(X0, F0);
      G0 = Gini;

      // Algorithm iterations
      real Tolerance_F_norm = this->m_Tolerance;
      real Tolerance_D_norm = this->m_Tolerance*this->m_Tolerance;
      this->m_Converged = false;
      for (
        this->m_Iterations = integer(1);
        this->m_Iterations < this->m_MaxIterations;
        ++this->m_Iterations
      )
      {
        // Calculate step
        this->step(F0, G0, D0);

        // Update history
        this->m_HistoryX.push_back(X0);
        this->m_HistoryF.push_back(F0);
        this->m_HistoryD.push_back(D0);

        // Check convergence
        F0_norm = F0.norm();
        D0_norm = D0.norm();
        if (F0_norm < Tolerance_F_norm || D0_norm < Tolerance_D_norm)
        {
          this->m_Converged = true;
          break;
        }

        // Update point
        X1 = X0 + D0;
        this->eval(X1, F1);

        // Update jacobian approximation
        DX1 = X1 - X0;
        DF1 = F1 - F0;
        this->update(
          X0, F0, D0, DX0, DF0, G0, // Old step data
          X1, F1, D1, DX1, DF1, G1  // New step data
        );

        // Print iteration data
        if (this->m_Verbose)
        {
          this->print(
            X0, F0, D0, DX0, DF0, G0, // Old step data
            X1, F1, D1, DX1, DF1, G1  // New step data
          );
        }

        // Update internal variables
        X0  = X1;
        F0  = F1;
        DX0 = DX1;
        DF0 = DF1;
        D0  = D1;
        G0  = G1;
      }

      // Convergence data
      Xout = X0;
      this->m_Residuals = F0_norm;
      return this->m_Converged;
    }

    //! Solve non-linear system of equations F(X)=0 with dumping
    bool
    solveDumped(
      vecT const & Xini, //!< Initialization point
      matT const & Gini, //!< Initialization jacobian approximation
      vecT       & Xout  //!< Output point
    )
    {
      // Setup internal variables
      this->reset();

      // Initialize variables
      matT G0, G1;
      vecT X0, F0, D0, DX0, DF0, X1, F1, D1, DX1, DF1;
      real F0_norm = real(0.0);
      real D0_norm = real(0.0);
      real F1_norm = real(0.0);
      real D1_norm = real(0.0);

      // Set initial iteration
      X0 = Xini;
      this->eval(X0, F0);
      G0 = Gini;

      // Algorithm iterations
      real Tolerance_F_norm = this->m_Tolerance;
      real Tolerance_D_norm = this->m_Tolerance*this->m_Tolerance;
      this->m_Converged = false;
      real tau = real(1.0);
      for (
        this->m_Iterations = integer(1);
        this->m_Iterations < this->m_MaxIterations;
        ++this->m_Iterations
      )
      {
        // Calculate step
        this->step(F0, G0, D0);

        // Update history
        this->m_HistoryX.push_back(X0);
        this->m_HistoryF.push_back(F0);
        this->m_HistoryD.push_back(D0);

        // Check convergence
        F0_norm = F0.norm();
        D0_norm = D0.norm();
        if (F0_norm < Tolerance_F_norm || D0_norm < Tolerance_D_norm)
        {
          this->m_Converged = true;
          break;
        }

        // Relax the iteration process
        tau = 1.0;
        for (
          this->m_Relaxations = integer(0);
          this->m_Relaxations < this->m_MaxRelaxations;
          ++this->m_Relaxations
        )
        {
          // Update point
          D1 = tau * D0;
          X1 = X0 + D1;
          this->eval(X1, F1);

          // Check relaxation
          F1_norm = F1.norm();
          D1_norm = D1.norm();
          if (
            F1_norm < F0_norm ||
            D1_norm < (real(1.0)-tau/real(2.0))*D0_norm
          )
            {break;}
          else
            {tau *= this->m_Alpha;}
        }

        // Update jacobian approximation
        DX1 = X1 - X0;
        DF1 = F1 - F0;
        this->update(
          X0, F0, D0, DX0, DF0, G0, // Old step data
          X1, F1, D1, DX1, DF1, G1  // New step data
        );

        // Print iteration data
        if (this->m_Verbose)
        {
          this->print(
            X0, F0, D0, DX0, DF0, G0, // Old step data
            X1, F1, D1, DX1, DF1, G1  // New step data
          );
        }

        // Update internal variables
        X0  = X1;
        F0  = F1;
        DX0 = DX1;
        DF0 = DF1;
        D0  = D1;
        G0  = G1;
      }

      // Convergence data
      Xout = X0;
      this->m_Residuals = F0_norm;
      return this->m_Converged;
    }

    //! Print iterations information
    void
    print(
      vecT const & X0,  //!< Old points
      vecT const & F0,  //!< Old function value
      vecT const & D0,  //!< Old step
      vecT const & DX0, //!< Old points difference
      vecT const & DF0, //!< Old function value difference
      matT const & G0,  //!< Old jacobian approximation
      vecT const & X1,  //!< New points
      vecT const & F1,  //!< New function value
      vecT const & D1,  //!< New step
      vecT const & DX1, //!< New points difference
      vecT const & DF1, //!< New function value difference
      matT const & G1   //!< New jacobian approximation
    )
    const
    {
      std::cout
        << "----" << std::endl
        << "IT = " << this->m_Iterations << std::endl
        << "RL = " << this->m_Relaxations << std::endl
        << "\t X0  = [ " << X0.transpose() << " ]'" << std::endl
        << "\t F0  = [ " << F0.transpose() << " ]'" << std::endl
        << "\t D0  = [ " << D0.transpose() << " ]'" << std::endl
        << "\t DX0 = [ " << DX0.transpose() << " ]'" << std::endl
        << "\t DF0 = [ " << DF0.transpose() << " ]'" << std::endl
        << "\t G0  = [ ";
        for (integer i = 0; i < G0.rows(); ++i)
          {std::cout << G0.row(i) << ";";}
      std::cout
        << "\b ]'" << std::endl
        << "\t X1  = [ " << X1.transpose()  << " ]'" << std::endl
        << "\t F1  = [ " << F1.transpose()  << " ]'" << std::endl
        << "\t D1  = [ " << D1.transpose()  << " ]'" << std::endl
        << "\t DX1 = [ " << DX1.transpose() << " ]'" << std::endl
        << "\t DF1 = [ " << DF1.transpose() << " ]'" << std::endl
        << "\t G1  = [ ";
        for (integer i = 0; i < G1.rows(); ++i)
          {std::cout << G1.row(i) << ";";}
      std::cout
        << "\b ]'" << std::endl
        << "----" << std::endl;
    }

    //! Jacobian approximation update rule
    virtual
    void
    update(
      vecT const & X0,  //!< Old points
      vecT const & F0,  //!< Old function value
      vecT const & D0,  //!< Old step
      vecT const & DX0, //!< Old points difference
      vecT const & DF0, //!< Old function value difference
      matT const & G0,  //!< Old jacobian approximation
      vecT const & X1,  //!< New points
      vecT const & F1,  //!< New function value
      vecT const & D1,  //!< New step
      vecT const & DX1, //!< New points difference
      vecT const & DF1, //!< New function value difference
      matT       & G1   //!< New jacobian approximation
    ) = 0;

    //! Calculate step
    virtual
    void
    step(
      vecT const & F, //!< Function
      matT const & G, //!< Jacobian approximation
      vecT       & D  //!< Step
    ) const = 0;

  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*\
   |   _____ _           _       _   _                       _ _
   |  | ____(_)_ __ ___ | | __ _| \ | | _____   ____ _ _ __ | (_)_ __  _ __   __ _
   |  |  _| | | '__/ _ \| |/ _` |  \| |/ _ \ \ / / _` | '_ \| | | '_ \| '_ \ / _` |
   |  | |___| | | | (_) | | (_| | |\  |  __/\ V / (_| | | | | | | | | | | | | (_| |
   |  |_____|_|_|  \___/|_|\__,_|_| \_|\___| \_/ \__,_|_| |_|_|_|_| |_|_| |_|\__,_|
   |
  \*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! QuasiNewton Eirola-Nevalinna method class container
  template <typename vecT, typename matT>
  class EirolaNevanlinna : public Solver<vecT, matT>
  {
  public:
    //! Return method name
    std::string
    method(void)
    const override
    {
      return "EirolaNevanlinna";
    }

  private:
    //! Jacobian approximation update rule
    void
    update(
      vecT const & X0,  //!< Old points
      vecT const & F0,  //!< Old function value
      vecT const & D0,  //!< Old step
      vecT const & DX0, //!< Old points difference
      vecT const & DF0, //!< Old function value difference
      matT const & G0,  //!< Old jacobian approximation
      vecT const & X1,  //!< New points
      vecT const & F1,  //!< New function value
      vecT const & D1,  //!< New step
      vecT const & DX1, //!< New points difference
      vecT const & DF1, //!< New function value difference
      matT       & G1   //!< New jacobian approximation
    )
    override
    {
      vecT FQ;
      vecT DQ(-G0 * F1);
      this->eval(X1 + DQ, FQ);
      G1 = G0 + (DQ - G0 * (FQ - F1)) * (DQ.transpose() * G0) / (DQ.transpose() * G0 * (FQ - F1));
    }

    //! Calculate step
    void
    step(
      vecT const & F, //!< Function
      matT const & G, //!< Jacobian approximation
      vecT       & D  //!< Step
    )
    const override
    {
      D = -G*F;
    }

  }; // class EirolaNevanlinna

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*\
   |   ____                      _
   |  | __ ) _ __ ___  _   _  __| | ___ _ __
   |  |  _ \| '__/ _ \| | | |/ _` |/ _ \ '_ \
   |  | |_) | | | (_) | |_| | (_| |  __/ | | |
   |  |____/|_|  \___/ \__, |\__,_|\___|_| |_|
   |                   |___/
  \*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! QuasiNewton Broyden's Ugly method class container
  template <typename vecT, typename matT>
  class BroydenUgly : public Solver<vecT, matT>
  {
  public:
    //! Return method name
    std::string
    method(void)
    const override
    {
      return "BroydenUgly";
    }

  private:
    //! Jacobian approximation update rule
    void
    update(
      vecT const & X0,  //!< Old points
      vecT const & F0,  //!< Old function value
      vecT const & D0,  //!< Old step
      vecT const & DX0, //!< Old points difference
      vecT const & DF0, //!< Old function value difference
      matT const & G0,  //!< Old jacobian approximation
      vecT const & X1,  //!< New points
      vecT const & F1,  //!< New function value
      vecT const & D1,  //!< New step
      vecT const & DX1, //!< New points difference
      vecT const & DF1, //!< New function value difference
      matT       & G1   //!< New jacobian approximation
    )
    override
    {
      // Broyden's Ugly method
      // G1 = G0 - (G0*DF1-DX1)/(C'*DX1)*C', where C = DX1;
      G1 = G0 - (G0*DF1-DF1)/(DX1.transpose()*DX1)*DX1.transpose();
    }

    //! Calculate step
    void
    step(
      vecT const & F, //!< Function
      matT const & G, //!< Jacobian approximation
      vecT       & D  //!< Step
    )
    const override
    {
      D = -G.inverse()*F;
    }

  }; // class BroydenUgly

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! QuasiNewton Broyden's Bad method class container
  template <typename vecT, typename matT>
  class BroydenBad : public Solver<vecT, matT>
  {
  public:
    //! Return method name
    std::string
    method(void)
    const override
    {
      return "BroydenBad";
    }

  private:
    //! Jacobian approximation update rule
    void
    update(
      vecT const & X0,  //!< Old points
      vecT const & F0,  //!< Old function value
      vecT const & D0,  //!< Old step
      vecT const & DX0, //!< Old points difference
      vecT const & DF0, //!< Old function value difference
      matT const & G0,  //!< Old jacobian approximation
      vecT const & X1,  //!< New points
      vecT const & F1,  //!< New function value
      vecT const & D1,  //!< New step
      vecT const & DX1, //!< New points difference
      vecT const & DF1, //!< New function value difference
      matT       & G1   //!< New jacobian approximation
    )
    override
    {
      // Broyden's Bad method
      // G1 = G0 - (G0*DF1-DX1)/(C'*DF1)*C', where C = DF1;
      G1 = G0 - (G0*DF1-DX1)/(DF1.transpose()*DF1)*DF1.transpose();
    }

    //! Calculate step
    void
    step(
      vecT const & F, //!< Function
      matT const & G, //!< Jacobian approximation
      vecT       & D  //!< Step
    )
    const override
    {
      D = -G*F;
    }

  }; // class BroydenBad

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! QuasiNewton Broyden's Good method class container
  template <typename vecT, typename matT>
  class BroydenGood : public Solver<vecT, matT>
  {
  public:
    //! Return method name
    std::string
    method(void)
    const override
    {
      return "BroydenGood";
    }

  private:
    //! Jacobian approximation update rule
    void
    update(
      vecT const & X0,  //!< Old points
      vecT const & F0,  //!< Old function value
      vecT const & D0,  //!< Old step
      vecT const & DX0, //!< Old points difference
      vecT const & DF0, //!< Old function value difference
      matT const & G0,  //!< Old jacobian approximation
      vecT const & X1,  //!< New points
      vecT const & F1,  //!< New function value
      vecT const & D1,  //!< New step
      vecT const & DX1, //!< New points difference
      vecT const & DF1, //!< New function value difference
      matT       & G1   //!< New jacobian approximation
    )
    override
    {
      // Broyden's Good method
      // G1 = G0 - (G0*DF1-DX1)/(C'*DF1)*C', where C = G0'*DX1;
      vecT C( G0.transpose()*DX1 );
      G1 = G0 - (G0*DF1-DX1)/(C.transpose()*DF1)*C.transpose();
    }

    //! Calculate step
    void
    step(
      vecT const & F, //!< Function
      matT const & G, //!< Jacobian approximation
      vecT       & D  //!< Step
    )
    const override
    {
      D = -G*F;
    }

  }; // class BroydenGood

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! QuasiNewton Broyden's Combined method class container
  template <typename vecT, typename matT>
  class BroydenCombined : public Solver<vecT, matT>
  {
  public:
    //! Return method name
    std::string
    method(void)
    const override
    {
      return "BroydenCombined";
    }

  private:
    //! Jacobian approximation update rule
    void
    update(
      vecT const & X0,  //!< Old points
      vecT const & F0,  //!< Old function value
      vecT const & D0,  //!< Old step
      vecT const & DX0, //!< Old points difference
      vecT const & DF0, //!< Old function value difference
      matT const & G0,  //!< Old jacobian approximation
      vecT const & X1,  //!< New points
      vecT const & F1,  //!< New function value
      vecT const & D1,  //!< New step
      vecT const & DX1, //!< New points difference
      vecT const & DF1, //!< New function value difference
      matT       & G1   //!< New jacobian approximation
    )
    override
    {
      vecT G0DF1( G0*DF1 );
      real DX1G0DF1 = DX1.transpose()*G0DF1;
      real DF1DF1   = DF1.transpose()*DF1;
      // Selection criteria
      // |DX1'*DX0|/|DX1'*G0*DF1| < |DF1'*DF0|/(DF1'*DF1)
      if (
        std::abs(DX1.transpose()*DX0)/std::abs(DX1G0DF1) <
        std::abs(DF1.transpose()*DF0)/DF1DF1 ||
        this->outIterations() < integer(2)
      )
      {
        // Broyden's Good method
        // G1 = G0 - (G0*DF1-DX1)/(C'*DF1)*C', where C = G0'*DX1;
        vecT C( G0.transpose()*DX1 );
        G1 = G0 - (G0DF1-DX1)/(C.transpose()*DF1)*C.transpose();
      }
      else
      {
        // Broyden's Bad method
        // G1 = G0 - (G0*DF1-DX1)/(C'*DF1)*C', where C = DF1;
        G1 = G0 - (G0DF1-DX1)/(DF1DF1)*DF1.transpose();
      }
    }

    //! Calculate step
    void
    step(
      vecT const & F, //!< Function
      matT const & G, //!< Jacobian approximation
      vecT       & D  //!< Step
    )
    const override
    {
      D = -G*F;
    }

  }; // class BroydenCombined

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*\
   |    ____                         _            _ _
   |   / ___|_ __ ___  ___ _ __  ___| |_ __ _  __| | |_
   |  | |  _| '__/ _ \/ _ \ '_ \/ __| __/ _` |/ _` | __|
   |  | |_| | | |  __/  __/ | | \__ \ || (_| | (_| | |_
   |   \____|_|  \___|\___|_| |_|___/\__\__,_|\__,_|\__|
   |
  \*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! QuasiNewton Greenstadt's 1st method class container
  template <typename vecT, typename matT>
  class Greenstadt1 : public Solver<vecT, matT>
  {
  public:
    //! Return method name
    std::string
    method(void)
    const override
    {
      return "Greenstadt1";
    }

  private:
    //! Jacobian approximation update rule
    void
    update(
      vecT const & X0,  //!< Old points
      vecT const & F0,  //!< Old function value
      vecT const & D0,  //!< Old step
      vecT const & DX0, //!< Old points difference
      vecT const & DF0, //!< Old function value difference
      matT const & G0,  //!< Old jacobian approximation
      vecT const & X1,  //!< New points
      vecT const & F1,  //!< New function value
      vecT const & D1,  //!< New step
      vecT const & DX1, //!< New points difference
      vecT const & DF1, //!< New function value difference
      matT       & G1   //!< New jacobian approximation
    )
    override
    {
      // Greenstadt's 1st method
      // G1 = G0 - (G0*DF1-DX1)/(C'*DF1)*C', where C = F1;
      G1 = G0 - (G0*DF1-DX1)/(F1.transpose()*DF1)*F1.transpose();
    }

    //! Calculate step
    void
    step(
      vecT const & F, //!< Function
      matT const & G, //!< Jacobian approximation
      vecT       & D  //!< Step
    )
    const override
    {
      D = -G*F;
    }

  }; // class Greenstadt1

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! QuasiNewton Greenstadt's 2nd method class container
  template <typename vecT, typename matT>
  class Greenstadt2 : public Solver<vecT, matT>
  {
  public:
    //! Return method name
    std::string
    method(void)
    const override
    {
      return "Greenstadt2";
    }

  private:
    //! Jacobian approximation update rule
    void
    update(
      vecT const & X0,  //!< Old points
      vecT const & F0,  //!< Old function value
      vecT const & D0,  //!< Old step
      vecT const & DX0, //!< Old points difference
      vecT const & DF0, //!< Old function value difference
      matT const & G0,  //!< Old jacobian approximation
      vecT const & X1,  //!< New points
      vecT const & F1,  //!< New function value
      vecT const & D1,  //!< New step
      vecT const & DX1, //!< New points difference
      vecT const & DF1, //!< New function value difference
      matT       & G1   //!< New jacobian approximation
    )
    override
    {
      // Greenstadt's 2nd method
      // G1 = G0 - (G0*DF1-DX1)/(C'*DF1)*C', where C  = G0'*G0*DF1;
      vecT C( G0.transpose()*G0*DF1 );
      G1 = G0 - (G0*DF1-DX1)/(C.transpose()*DF1)*C.transpose();
    }

    //! Calculate step
    void
    step(
      vecT const & F, //!< Function
      matT const & G, //!< Jacobian approximation
      vecT       & D  //!< Step
    )
    const override
    {
      D = -G*F;
    }

  }; // class Greenstadt2

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

} // namespace QuasiNewton

#endif

///
/// eof: Solver.hxx
///
