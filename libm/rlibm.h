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


double rlibm_log10f(float);
double rlibm_log2f(float);
double rlibm_logf(float);
double rlibm_exp10f(float);
double rlibm_exp2f(float);
double rlibm_expf(float);
double rlibm_coshf(float);
double rlibm_sinhf(float);
double rlibm_cospif(float);
double rlibm_sinpif(float);
double rlibm_sinf(float);
double rlibm_cosf(float);
double rlibm_tanf(float);
double rlibm_atanf(float);
double rlibm_asinf(float);
double rlibm_acosf(float);

#endif
