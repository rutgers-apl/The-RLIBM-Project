#include "mpfr.h"
#include "rlibm.h"
#include <stdio.h>
#include <stdlib.h>
#include<assert.h>

#include "rlibm_log.h"

int main(int argc, char** argv){

  mpfr_t mval;
  mpfr_t mval2;

  mpfr_init2(mval, 200);
  mpfr_init2(mval2, 200);

  double lut_Log2F[128];
  double lut_OneByF[128];

  for (unsigned int i  = 0; i< 128; i++){
    float_x F = {.x =i<<16};
    F.x = F.x | 0x3F800000;

    int status = mpfr_set_d(mval, (double)F.f, MPFR_RNDN);
    if(status != 0){
      printf("Something went wrong while setting float to mpfr\n");
      exit(0);
    }

    status = mpfr_log2(mval, mval, MPFR_RNDN);
    double Log2F = mpfr_get_d(mval, MPFR_RNDN);
    status = mpfr_set_d(mval, (double)F.f, MPFR_RNDN);
    
    if(status != 0){
      printf("Something went wrong while setting float to mpfr\n");
      exit(0);
    }

    status = mpfr_set_d(mval2, (double)1.0F, MPFR_RNDN);
    if(status != 0){
      printf("Something went wrong while setting float to mpfr\n");
      exit(0);
    }

    status = mpfr_div(mval, mval2, mval, MPFR_RNDN);

    double OnebyLog2F = mpfr_get_d(mval, MPFR_RNDN);
    //    printf("F=%a, Log2F = %a, OnebyLog2F=%a\n", F.f,Log2F, OnebyLog2F);   

    lut_Log2F[i] = Log2F;
    lut_OneByF[i] = OnebyLog2F;

    double_x result1 = {.d = rlibm_log2F[i]};
    double_x result2 = {.d = Log2F};
    assert(result1.x == result2.x);

    result1.d = rlibm_OneByF[i];
    result2.d = OnebyLog2F;
    assert(result1.x == result2.x); 
  }
  //  printf("tables for Log2F:\n");
  printf("double rlibm_log2F[128]={\n");
  for(int i = 0; i<128; i++){
    if(i != 127){
      printf("%a, \n", lut_Log2F[i]);
    }
    else{
      printf("%a \n}\n", lut_Log2F[i]);
    }  
  }

  //  printf("table for OneByF:\n");
  printf("double rlibm_OneByF[128]={\n");
  for(int i = 0; i<128; i++){
    if(i != 127){
      printf("%a, \n", lut_OneByF[i]);
    }
    else{
      printf("%a \n}\n", lut_OneByF[i]);
    }  
  }
  mpfr_clear(mval);
  mpfr_clear(mval2);
  return 0;
}
