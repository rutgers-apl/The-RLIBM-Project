# The RLIBM Project

RLIBM is a collection of elementary functions that produces correctly
rounded results for multiple representations for multiple rounding
modes with a single implementation. The RLIBM project makes a case for
approximating the correctly rounded result of an elementary function
rather than the real value of an elementary function. When we
approximate the correctly rounded result, there is an interval of real
values around the correctly rounded result such that producing a real
value in this interval rounds to the correct result. This interval is
the freedom that the polynomial approximation has for an input, which
is larger than the freedom with prior approaches. Hence, the RLIBM
approach has more margin to generate correct, yet, efficient
polynomials.

If you are interested, see more details about the [RLIBM
project](https://people.cs.rutgers.edu/~sn349/rlibm/)

## Building the latest RLIBM implementations

The latest RLIBM functions that have been tested to produce correct results
for all representations from 10-bits to 32-bits with all the five
rounding modes (FE_NEAREST, FE_UPWARD, FE_DOWNWARD, FE_TOWARDZERO) is available in the libm folder.

The RLIBM library provides correctly rounded implementations for the following functions:

```
double rlibm_log10(float);
double rlibm_log2(float);
double rlibm_log(float);
double rlibm_exp10(float);
double rlibm_exp2(float);
double rlibm_exp(float);
double rlibm_cosh(float);
double rlibm_sinh(float);
double rlibm_cospi(float);
double rlibm_sinpi(float);
double rlibm_sinf(float);
double rlibm_cosf(float);
double rlibm_tanf(float);
```



### Buidling the latest RLIBM functions

```
cd libm
make 
```

### Using the RLIBM functions

The RLIBM project implements all elementary functions in double
precision with round-to-nearest-ties-to-even mode. It returns a double
value that when rounded to any target representation or rounding mode
produces correct rounded results. For example, the prototype of log2
is as follows.

```
double rlibm_log2(float);
```

To produce correctly rounded results for the
round-to-nearest-ties-to-even mode, the recommended usage in the
application is as follows.

```
float x = .. ;
float y = (float) rlibm_log2(x);

```


To produce correctly rounded results for other rounding modes such as
round-towards-positive-infinity, the recommended usage in the
application is as follows.


```
float x = .. ;
double temp  = rlibm_log2(x);
fesetround(FE_UPWARD);
float y = (float)(temp);
fesetround(FE_TONEAREST);
```

Link the program with rlibm.a in libm/ folder.

# Generating Correctly Rounded Implementations

The RLIBM functions have been generated using our prior prototypes:
[RLIBM-ALL](https://github.com/rutgers-apl/rlibm-all) and
[RLIBM-PROG](https://github.com/santoshn/rlibm-prog)

Specifically, we use the polynomial generation algorithm from
[RLIBM-PROG, PLDI
2022](https://people.cs.rutgers.edu/~sn349/papers/rlibm-prog-pldi-2022.pdf)
and we approximate the correctly rounded round-to-odd result, which
allows us to produce correctly rounded results for multiple
representations with a single polynomial approximation as with
[RLIBM-ALL, POPL
2022](https://people.cs.rutgers.edu/~sn349/papers/rlibmall-popl-2022.pdf).
