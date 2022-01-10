# The-RLIBM-Project
[Work in Progress]
See details about the RLIBM project here at https://people.cs.rutgers.edu/~sn349/rlibm/

# Building RLIBM-ALL functions

The RLIBM-ALL functions are in libm/rlibm-prog/libm folder.

1. Building RLIBM-ALL functions

```
cd libm/rlibm-prog/libm
make floatrnolibm
```

2. Link the test harness or the functions with floatrnolibm.a in libm/rlibm-prog/libm folder

3. Function prototype of the RLIBM-ALL functions for producing correctly rounded results for all rounding modes for a 32-bit float:

```
float rlibm_all_fast_log2(float);
float rlibm_all_fast_log10(float);
float rlibm_all_fast_log(float);
float rlibm_all_fast_exp2(float);
float rlibm_all_fast_exp10(float);
float rlibm_all_fast_exp(float);
float rlibm_all_fast_sinpi(float);
float rlibm_all_fast_cospi(float);
float rlibm_all_fast_sinh(float);
float rlibm_all_fast_cosh(float);
```

4.  To generate the correctly rounded result of a 32-bit float for a specific rounding
   mode such as FE_TONEAREST, FE_UPWARD, FE_DOWNWARD, or
   FE_TOWARDZERO, use fesetround function as follows:

```
fesetround(FE_TOWARDZERO);
float result = rlibm_all_fast_log(x);
```

