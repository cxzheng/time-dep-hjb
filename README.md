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

#### Parallel Computation
The code exploit multi-processor CPU to enable parallel compution (using __OPENMP__). The numerical solvers can be accelerated by parallel computation especially when the grid resolution is high. To enable parallel computation, turn on the flag `USE_OPENMP` in __cmake__.

### Customizing the PDEs
The code solves the time-dependent Hamilton-Jacobi-Bellman (HJB) PDEs which have the form

![$$u_t + f(x,t)|\nabla u| = g(x,t).$$](https://raw.githubusercontent.com/cxzheng/time-dep-hjb/master/images/jqjt6g4.png)

This code includes three different numerical solvers for this type of PDE, namely the explicit, implicit, and hybrid methods. The main C++ files for these solvers are `ExplicitTDHJ2D.cpp`, `ImplicitTDHJ2D.cpp`, and `MixedTDHJ2D.cpp`, respectively. 
These solvers are implemented as C++ templates, with both _f_ and _g_ functions as well as the boundary conditions specified 
as the template parameters.

Both the _f_ and _g_ functions are defined as C++ functionals, (i.e., a C++ struct with the operator `()`). For example, a _f_ function is defined as
```C++
struct func_F
{
    inline double operator() (const vector2d& pos, double t) const
    {   return ... }
};
```
And the _g_ function is defined in a similar way (e.g., see `src/HJ2DFunc.h`),
```C++
struct func_G
{
    inline double operator() (const vector2d& pos, double t) const
    {   return ... }
};
```
The boundary condition class specifies the number of boundary nodes, the positions of the boundary nodes, and their values. Please refoer to the code in `src/BoundaryCond.h` for a few examples.

Once _f_, _g_ and the boundary condition class are defined, you can solve the PDE by defining a test case in the solver’s main code (e.g. `ImplicitTDHJ2D.cpp`),
```C++
#ifdef USE_TEST_XX
typedef SimpleBoundaryCond2D                      TBC;
typedef CachedImpLattice2D<func_F, func_G, TBC>   TLattice;
```
Here `USE_TEST_XX` is the test case number. You can now enable this test case by define a flag `USE_TEST_XX` in `src/config.h` or add it in the compile flags.

#### Example
Turn on the testcase 07 (which is the experiment 4 in the paper) by providing a `config.h' in `src` directory
```C++
#ifndef CONFIG_INC
#   define CONFIG_INC
#define  USE_TEST_07
#endif
```
After compilation, run 
```
src/ExplicitTDHJ2D -X 128 -Y 128 -T 4 -d 0.001104 -o u128_0.001104_4.npy
```
to solve the PDE with a 128x128 grid and a timestep size of `0.001104`, using explicit method. The simulation runs backward from `T=4` down to `T=0`.

### Interpreting the Results
In the above example, the final data __u(x,0)__ is stored in the file `u128_0.001104_4.npy`. This files stores the 2D array of data at time slice __t=0__, and the file format is Numpy's [NPY format](https://docs.scipy.org/doc/numpy-dev/neps/npy-format.html), so it can be easily read and visualized by python. For example, one can visualize the result using
```
../python/npy_plot.py u128_0.001104_4.npy
```
In addition, the folder `test_scripts` includes two Perl scripts as an example of testing the testcase 07 using explicit, implicit, and hybrid methods with various resolutions. In particular, `test_scripts/test07.pl` runs the numerical solves, and `test_scripts/test07_err.pl` evaluates the accuracy of different solvers by comparing with a high-resolution (3072x3072) result (using explicit method).
