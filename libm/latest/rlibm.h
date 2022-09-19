#ifndef __RLIBM_LATEST_LIBM_H__
#define __RLIBM_LATEST_LIBM_H__

#include <stdint.h>

typedef union {
  double d;
  uint64_t x;
} double_x;

typedef union {
  float f;
  uint32_t x;
} float_x;


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

#endif
