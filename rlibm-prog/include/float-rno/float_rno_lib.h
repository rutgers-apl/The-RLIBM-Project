#ifndef __FLOAT_ROUND_TO_ODD_LIB_H__
#define __FLOAT_ROUND_TO_ODD_LIB_H__

#include "common.h"
#include <stdint.h>


double rlibm_rno_log10(float);
double rlibm_rno_log10_piecewise(float);
double rlibm_rno_log10_piecewise_v2(float);
double rlibm_rno_log2(float);
double rlibm_rno_log2_piecewise(float);
double rlibm_rno_log(float);
double rlibm_rno_log_piecewise(float);
double rlibm_rno_exp10(float);
double rlibm_rno_exp10_piecewise(float);
double rlibm_rno_exp2(float);
double rlibm_rno_exp2_piecewise(float);
double rlibm_rno_exp(float);
double rlibm_rno_exp_piecewise(float);
double rlibm_rno_cosh(float);
double rlibm_rno_sinh(float);
double rlibm_rno_cospi(float);
double rlibm_rno_sinpi(float);

#if 0
float rlibm_rno_exp2(float);

float rlibm_fast_exp(float);
float rlibm_fast_log10(float);
float rlibm_fast_exp10(float);
float rlibm_fast_exp10_v2(float);
float rlibm_fast_log(float);
float rlibm_fast_cosh(float);
float rlibm_fast_cosh_v2(float);
float rlibm_fast_cosh_v3(float);
float rlibm_fast_cosh_v4(float);
float rlibm_fast_sinh(float);
float rlibm_fast_cospi(float);
float rlibm_fast_sinpi(float);

// rlibm32 headers
float rlibm_log(float);
float rlibm_log2(float);
float rlibm_log10(float);
float rlibm_exp(float);
float rlibm_exp2(float);
float rlibm_exp10(float);
float rlibm_sinh(float);
float rlibm_cosh(float);
float rlibm_sinpi(float);
float rlibm_cospi(float);
float rlibm_log2_8(float);
float rlibm_log10_8(float);
#endif

#endif
