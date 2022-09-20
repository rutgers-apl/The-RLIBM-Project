#include <math.h>
#include "float-rno/float_rno_lib.h"
#include <fenv.h>

float rlibm_all_fast_log2(float x){

  int prev_round_mode = fegetround();
  fesetround(FE_TONEAREST);
  double val = rlibm_rno_log2(x);
  fesetround(prev_round_mode);
  float result = (float) val;
  return result;
}

float rlibm_all_fast_log(float x){

  int prev_round_mode = fegetround();
  fesetround(FE_TONEAREST);
  double val = rlibm_rno_log_piecewise(x);
  fesetround(prev_round_mode);
  float result = (float) val;
  return result;
}

float rlibm_all_fast_log10(float x){

  int prev_round_mode = fegetround();
  fesetround(FE_TONEAREST);
  double val = rlibm_rno_log10_piecewise(x);
  fesetround(prev_round_mode);
  float result = (float) val;
  return result;
}

float rlibm_all_fast_exp2(float x){

  int prev_round_mode = fegetround();
  fesetround(FE_TONEAREST);
  double val = rlibm_rno_exp2(x);
  fesetround(prev_round_mode);
  float result = (float) val;
  return result;
}

float rlibm_all_fast_exp(float x){

  int prev_round_mode = fegetround();
  fesetround(FE_TONEAREST);
  double val = rlibm_rno_exp_piecewise(x);
  fesetround(prev_round_mode);
  float result = (float) val;
  return result;  
}

float rlibm_all_fast_exp10(float x){

  int prev_round_mode = fegetround();
  fesetround(FE_TONEAREST);
  double val = rlibm_rno_exp10_piecewise(x);
  fesetround(prev_round_mode);
  float result = (float) val;
  return result;  
}

float rlibm_all_fast_sinpi(float x){

  int prev_round_mode = fegetround();
  fesetround(FE_TONEAREST);
  double val = rlibm_rno_sinpi(x);
  fesetround(prev_round_mode);
  float result = (float) val;
  return result;  
}

float rlibm_all_fast_cospi(float x){

  int prev_round_mode = fegetround();
  fesetround(FE_TONEAREST);
  double val = rlibm_rno_cospi(x);
  fesetround(prev_round_mode);
  float result = (float) val;
  return result;  
}

float rlibm_all_fast_sinh(float x){

  int prev_round_mode = fegetround();
  fesetround(FE_TONEAREST);
  double val = rlibm_rno_sinh(x);
  fesetround(prev_round_mode);
  float result = (float) val;
  return result;  
}

float rlibm_all_fast_cosh(float x){

  int prev_round_mode = fegetround();
  fesetround(FE_TONEAREST);
  double val = rlibm_rno_cosh(x);
  fesetround(prev_round_mode);
  float result = (float) val;
  return result;  
}




