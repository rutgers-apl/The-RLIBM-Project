# Generating Polynomials

This folder contains the polynomial generator that generates
polynomials with specific terms given the interval file and the
configuration files. The polynomial generation needs SOPLEX-4.0.1
installed and the following environment variables set.


```
export SOPLEX_INCLUDE=<src directory of the SOPLEX-4.0.1>
export SOPLEX_LIB=<full path to libsoplex.a>
```
SOPLEX can be downloaded from [here](https://github.com/scipopt/soplex/releases/tag/release-401).


## Providing the polynomial evaluation

The polynomial generator generates polynomial assuming the specific
kind of polynomial evaluation (Horner's method, Estrin's method, or
Estrin+ FMA).

Provide the polynomial evaluation as function that returns a double
value given a input to polynomial evaluation and a polynomial. The
prototype of the function is 

```
double rlibm_poly_evaluation(double, polynomial*);
```


## Building the Polynomial Generator

To build the polynomial generator, execute the following command.

```
make
```

## Polynomial configuration file

The polynomial generator generates polynomials with the specified
features (piecewise, specific degrees, and terms). The first line of
the configuration file is the number of pieces. The second line
provides the number of terms and the degree of each term of the
polynomial.

```
1
5 1 2 3 4 5
```

## Running the polynomial generator

The polynomial generator needs the interval file and the configuration
file. 

```
./polygen <interval file> <polynomial configuration file>
```


