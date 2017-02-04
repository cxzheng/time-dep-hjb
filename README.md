# time-dep-hjb
[A fast implicit method for time-dependent Hamilton-Jacobi PDEs](http://www.cs.columbia.edu/~cxz/TimeDepHJB/)

## Abstract
We present a new efficient computational approach for time-dependent first-order Hamilton-Jacobi-Bellman PDEs. Since our method is based on a time-implicit Eulerian discretization, the numerical scheme is unconditionally stable, but discretized equations for each time-slice are coupled and non-linear. We show that the same system can be re-interpreted as a discretization of a static Hamilton-Jacobi-Bellman PDE on the same physical domain. The latter was shown to be "causal", making fast (non-iterative) methods applicable. The implicit discretization results in higher computational cost per time slice compared to the explicit time marching. However, the explicit discretization is subject to a CFL-stability condition, and the implicit approach becomes significantly more efficient whenever the accuracy demands on the time-step are less restrictive than the stability. We also present a hybrid method, which aims to combine the advantages of both the explicit and implicit discretizations. We demonstrate the efficiency of our approach using several 2D examples in optimal control of isotropic fixed-horizon processes.

## About the Source Code
This is a C++ implementation of the paper [A fast implicit method for time-dependent Hamilton-Jacobi PDEs](http://www.cs.columbia.edu/~cxz/TimeDepHJB/) by Vladimirsky and Zheng.

### Code Compilation
The code compilation has been tested in Linux using both __gcc__ and Intel’s C++ compiler (__icpc__). If you have CMake installed, you can compile it using the following steps:

* Create a buld directory and go into it
```
mkdir gcc-build && cd gcc-build
```
* Run __cmake__ to configure the Makefile
```
cmake ..
```
* Compile the code
```
make
```
After compiling the code, you can run the three executables located in gcc-build/src folder. They are __ExplicitTDHJ2D__, __ImplicitTDHJ2D__ and __MixedTDHJ2D__ respectively corresponding to the explicit, implicit and hybrid methods described in the paper. Different command line options are possible. See their details by providing the `-h` option (e.g., run `ImplicitTDHJ2D -h`).

#### Required Library
Compiling this source code requires [Boost](http://www.boost.org/) library.

### Customizing the PDEs
The code solves the time-dependent Hamilton-Jacobi-Bellman (HJB) PDEs which have the form

![$$u_t + f(x,t)|\nabla u| = g(x,t).$$](https://raw.githubusercontent.com/cxzheng/time-dep-hjb/master/images/jqjt6g4.png)

This code includes three different numerical solvers for this PDE.
