# Project for the PhD course of "Scientific Programming"

Course name: *Scientific Programming*
Professor name: *Enrico Bertolazzi*
PhD Student: *Davide Stocco*
Academic year: *2022/2023*

## Introduction

This repository contains a small collection of Quasi-Newton solvers. In particular, the methods that are implemented are the following:

  - Greenstadt’s 1st Method;
  - Greenstadt’s 2nd Method;
  - Eirola-Nevanlinna Method
  - Broyden’s Bad Method;
  - Broyden’s Good Method;
  - Broyden’s Combined Method;

and their respective dumped versions (indicated with the prefix “D–XXM”). We will compare the performance of the algorithms in terms of function evaluations, iterations, convergence speed, and success ratio on a set of tests.

## Implementation

All methods implemented in this small library have the same frame. In particular, an abstract class ``Solver`` is implemented. This class contains all the evaluation counters, the damping parameter, and the ``solve()`` and ``solveDumped()`` methods. On the other hand, the methods ``method()``, ``step()``, and ``update()`` are declared as virtual and will later be overridden in the specialised classes that will implement the above solvers. The ``Solver`` class also has all the setters and getters to trick it.

### Implementation

We will make use of the following quantities:

$$
\begin{align}
F_{k}        &= F_{k}(x_{k})    \\
\Delta x_{k} &= x_{k+1} - x_{k} \\
\Delta F_{k} &= F_{k+1} - F_{k} \\
\end{align}
$$

#### Greenstadt’s 1st Method (G1M)

For the *Greenstadt’s 1st Method*, the step for the inverse approximate Jacobian matrix update is:

$$
\begin{align}
G_{k+1}^{-1} &= G_{k}^{-1} - {(G_{k}^{-1} {\Delta F_{k}} - {\Delta x_{k}}) C_{k}^T} \over C^T {\Delta F_{k}} \\
C_{k}        &= F_{k}
\end{align}
$$

#### Greenstadt’s 2nd Method (G2M)

For the *Greenstadt’s 2nd Method*, the step for the inverse approximate Jacobian matrix update is:

$$
\begin{align}
G_{k+1}^{-1} &= G_{k}^{-1} - {(G_{k}^{-1} {\Delta F_{k}} - {\Delta x_{k}}) C_{k}^T} \over C^T {\Delta F_{k}} \\
C_{k}        &= G_{k}^{-1\,T} G_{k}^{-1} \Delta F_{k}
\end{align}
$$

#### Broyden’s Ugly Method (BUM)

For the *Broyden’s Ugly Method*, the step for the approximate Jacobian matrix update is:

$$
\begin{align}
G_{k+1} &= G_{k} - {(G_{k} {\Delta F_{k}} - {\Delta x_{k}}) C_{k}^T} \over C^T {\Delta F_{k}} \\
C_{k}   &= G_{k}^{T} G_{k} \Delta x_{k}
\end{align}
$$

#### Broyden’s Bad Method (BBM)

For the *Broyden’s Bad Method*, the step for the inverse approximate Jacobian matrix update is:

$$
\begin{align}
G_{k+1}^{-1} &= G_{k}^{-1} - {(G_{k}^{-1} {\Delta F_{k}} - {\Delta x_{k}}) C_{k}^T} \over C^T {\Delta F_{k}} \\
C_{k}        &= \Delta F_{k}
\end{align}
$$

#### Broyden’s Good Method (BGM)

For the *Broyden’s Good Method*, the step update for the inverse approximate Jacobian matrix is:

$$
\begin{align}
G_{k+1}^{-1} &= G_{k}^{-1} - {(G_{k}^{-1} {\Delta F_{k}} - {\Delta x_{k}}) C_{k}^T} \over C^T {\Delta F_{k}} \\
C_{k}        &= G_{k}^{-1\,T} \Delta x_{k}
\end{align}
$$

#### Broyden’s Combined Method (BCM)

For the *Broyden’s Combined Method*, the step for the inverse approximate Jacobian matrix update is:

$$
\begin{align}
G_{k+1}^{-1} &= G_{k}^{-1} - {{(G_{k}^{-1} {\Delta F_{k}} - {\Delta x_{k}}) C_{k}^T} \over {C^T {\Delta F_{k}}}} \\
C_{k} &=
\begin{cases}
G_{k}^{-1\,T} \Delta x_{k} & \mathrm{(BGM)} \quad \displaystyle{{|x_{k}^T x_{k-1}|} \over {|x_{k}^T G_{k}^{-1} x_{k}|}} < {|{F_{k}^T F_{k-1}|} \over {F_{k}^T F_{k}}} \\
\Delta F_{k}               & \mathrm{(BBM)} \quad \mathrm{otherwise}
\end{cases}
\end{align}
$$

#### Eirola-Nevanlinna Method (ENM)

For the *Eirola-Nevanlinna Method*, the step for the inverse approximate Jacobian matrix update is:

$$
\begin{align}
G_{k+1}^{-1} &= G_{k}^{-1} + {(P_k - G_{k}^{-1} Q_{k}) P_k^T G_{k}^{-1} \over P_k^T G_{k}^{-1} Q_{k}} \\
P_k   &= -G_{k-1} F_{k} \\
Q_{k} &= F(x_{k} + P_k) - F_{k}
\end{align}
$$

## References

  - Broyden
  - Broyden
  -

