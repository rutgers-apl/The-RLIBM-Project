#ifndef __FLOAT_PROGRESSIVE_LIB_H__
#define __FLOAT_PROGRESSIVE_LIB_H__

#include "common.h"
#include <stdint.h>

double rlibm_prog_rno_exp2(float);
double rlibm_prog_bf16_exp2(float);
double rlibm_prog_tf32_exp2(float);

double rlibm_prog_rno_exp(float);
double rlibm_prog_bf16_exp(float);
double rlibm_prog_tf32_exp(float);

double rlibm_prog_rno_cosh(float);
double rlibm_prog_bf16_cosh(float);
double rlibm_prog_tf32_cosh(float);

double rlibm_prog_rno_cospi(float);
double rlibm_prog_bf16_cospi(float);
double rlibm_prog_tf32_cospi(float);

double rlibm_prog_rno_sinh(float);
double rlibm_prog_bf16_sinh(float);
double rlibm_prog_tf32_sinh(float);

double rlibm_prog_rno_sinpi(float);
double rlibm_prog_bf16_sinpi(float);
double rlibm_prog_tf32_sinpi(float);

double rlibm_prog_rno_exp10(float);
double rlibm_prog_bf16_exp10(float);
double rlibm_prog_tf32_exp10(float);

double rlibm_prog_rno_log2(float);
double rlibm_prog_bf16_log2(float);
double rlibm_prog_tf32_log2(float);

double rlibm_prog_rno_log10(float);
double rlibm_prog_bf16_log10(float);
double rlibm_prog_tf32_log10(float);

double rlibm_prog_rno_log(float);
double rlibm_prog_bf16_log(float);
double rlibm_prog_tf32_log(float);

#endif
