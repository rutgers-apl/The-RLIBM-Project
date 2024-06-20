# Maximum Consensus Solver with FP Solutions for RLIBM Constraints

Given a set of infeasible (or feasible) low-dimensional linear
constraints, the maximum consensus solver produces a FP solution in
double precision that satisifies a maximum number of constraints. It
computes the convex hull in 2-D dimensions to identify a superset of
violated constraints, splits the system into a feasible subset and a
set of violated constraints. The implementation of the Clarkson's
method from the RLIBM project is used to identify the basis of the
feasible set. Subsequently, a new linear program that satisifes a
maximum number of constraints in the violated set while satsifying all
the constraints in the feasible set. See more details about this
appproach [here](https://people.cs.rutgers.edu/~sn349/papers/maxfs-pldi-2024.pdf). 

## Dependencies

To build the MCS solver, we need the following packages to be installed.

* GMP - GNU Multiple Precision Arithmetic Library
* MPFR - The GNU MPFR library
* CGAL

You can install these on a Debian Linux distribution with the following command

``` 
sudo apt-get install libgmp-dev libmpfr-dev libcgal-dev     
```

You also need SOPLEX-4.0.1 installed and the following environment
variables set.


```
export SOPLEX_INCLUDE=<src directory of the SOPLEX-4.0.1>
export SOPLEX_LIB=<full path to libsoplex.a>
```
SOPLEX can be downloaded from [here](https://github.com/scipopt/soplex/releases/tag/release-401).

## Building the MCS-Solver



## Computing the Convex Hull



### Building the Convex hull and Max Consensus Solver


