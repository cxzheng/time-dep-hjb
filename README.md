# time-dep-hjb
[A fast implicit method for time-dependent Hamilton-Jacobi PDEs](http://www.cs.columbia.edu/~cxz/TimeDepHJB/)

## Abstract
We present a new efficient computational approach for time-dependent first-order Hamilton-Jacobi-Bellman PDEs. Since our method is based on a time-implicit Eulerian discretization, the numerical scheme is unconditionally stable, but discretized equations for each time-slice are coupled and non-linear. We show that the same system can be re-interpreted as a discretization of a static Hamilton-Jacobi-Bellman PDE on the same physical domain. The latter was shown to be "causal", making fast (non-iterative) methods applicable. The implicit discretization results in higher computational cost per time slice compared to the explicit time marching. However, the explicit discretization is subject to a CFL-stability condition, and the implicit approach becomes significantly more efficient whenever the accuracy demands on the time-step are less restrictive than the stability. We also present a hybrid method, which aims to combine the advantages of both the explicit and implicit discretizations. We demonstrate the efficiency of our approach using several 2D examples in optimal control of isotropic fixed-horizon processes.

## Source Code README
This is a C++ implementation of the paper 
