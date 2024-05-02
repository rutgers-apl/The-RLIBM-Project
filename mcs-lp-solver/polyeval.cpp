#include "maxfs_common.h"

extern my_function function_to_process;


void set_function_process(char** argv){
  if(strcmp(argv[5], "log") == 0){
    function_to_process = LOG;
  }
  if(strcmp(argv[5], "log2") == 0){
    function_to_process = LOG2;
  }
  
  if(strcmp(argv[5], "log10") == 0){
    function_to_process = LOG10;
  }
  
  if(strcmp(argv[5], "exp2") == 0){
    function_to_process = EXP2;
  }
  
  if(strcmp(argv[5], "exp") == 0){
    function_to_process = EXP;
  }
  
  if(strcmp(argv[5], "exp10") == 0){
    function_to_process = EXP10;
  }

  if(strcmp(argv[5], "sinh_sinh") == 0){
    function_to_process = SINH_SINH;
  }

  if(strcmp(argv[5], "sinh_cosh") == 0){
    function_to_process = SINH_COSH;
  }

  if(strcmp(argv[5], "cosh_sinh") == 0){
    function_to_process = COSH_SINH;
  }

  if(strcmp(argv[5], "cosh_cosh") == 0){
    function_to_process = COSH_COSH;
  }

    if(strcmp(argv[5], "sinpi_sinpi") == 0){
    function_to_process = SINPI_SINPI;
  }

  if(strcmp(argv[5], "sinpi_cospi") == 0){
    function_to_process = SINPI_COSPI;
  }

  if(strcmp(argv[5], "cospi_sinpi") == 0){
    function_to_process = COSPI_SINPI;
  }

  if(strcmp(argv[5], "cospi_cospi") == 0){
    function_to_process = COSPI_COSPI;
  }

  if(strcmp(argv[5], "sin_sin") == 0){
    function_to_process = SIN_SIN;
  }

  if(strcmp(argv[5], "sin_cos") == 0){
    function_to_process = SIN_COS;
  }

  if(strcmp(argv[5], "cos_sin") == 0){
    function_to_process = COS_SIN;
  }

  if(strcmp(argv[5], "cos_cos") == 0){
    function_to_process = COS_COS;
  }



  assert(function_to_process != NONE);
}

double rlibm_poly_evaluation(double x, polynomial* poly){

  if(function_to_process == SINH_COSH || function_to_process == COSH_COSH || function_to_process == SINPI_COSPI || function_to_process == COSPI_COSPI || function_to_process == COS_COS || function_to_process == SIN_COS){
    assert(poly->termsize == 3);

    double C0 = poly->coeffs[0];
    double C2 = poly->coeffs[1];
    double C4 = poly->coeffs[2];

    double xsquare = x * x;

    double temp1 = fma(xsquare, C4, C2);
    return fma(xsquare, temp1, C0);    
  }
  
  if(function_to_process== SINH_SINH || function_to_process == COSH_SINH || function_to_process == SINPI_SINPI || function_to_process == COSPI_SINPI || function_to_process == SIN_SIN || function_to_process == COS_SIN){

    assert(poly->termsize == 3);

    double C1 = poly->coeffs[0];
    double C3 = poly->coeffs[1];
    double C5 = poly->coeffs[2];

    double xsquare = x * x;

    double temp1 = fma(xsquare, C5, C3);
    double temp2 = fma(xsquare, temp1, C1);
    return x * temp2;
  }


  if(function_to_process == LOG || function_to_process == LOG2 || function_to_process == LOG10){  

    assert(poly->termsize == 6);
    double  C1 = poly->coeffs[0];
    double  C2 = poly->coeffs[1];
    double  C3 = poly->coeffs[2];
    double  C4 = poly->coeffs[3];
    double  C5 = poly->coeffs[4];
    double  C6 = poly->coeffs[5];
    
    double xsquare = x * x;
    double xcube = x * xsquare;
    double temp1 = fma(C2, x, C1);
    //  double temp1 = C2 * x +  C1;
    double temp2 = fma(C3, xsquare, temp1);
    
    //  double temp2 = C3 * xsquare + temp1;
    
    double temp3 = fma(C5, x, C4);
    
    //  double temp3 = C5 * x + C4;
    
    double temp4 = fma(C6, xsquare, temp3);
    //  double temp4 = C6 * xsquare + temp3;
    
    double temp5 = fma(temp4, xcube, temp2);
    //  double temp5 = temp4 *xcube + temp2;
    
    return x * temp5;
  }

  if(function_to_process == EXP || function_to_process == EXP2 || function_to_process == EXP10){

    assert(poly->termsize == 6);
    double C0 = poly->coeffs[0];
    double C1 = poly->coeffs[1];
    double C2 = poly->coeffs[2];
    double C3 = poly->coeffs[3];
    double C4 = poly->coeffs[4];
    double C5 = poly->coeffs[5];
    
    double xsquare = x * x;
    double temp1 = fma (x, C1, C0);
    double temp2 = fma (x, C5, C4);
    double temp3 = fma (x, C3, C2);

    double temp4 = fma(xsquare, temp2, temp3);
    double temp5 = fma(xsquare, temp4, temp1);

    return temp5;
   
  }
  assert(0);
  return 0.0;
#if 0
  double ret_val = 0.0;

  // simulated Horner's method
  for(int i = poly->termsize-1; i> 0; i--){
    ret_val = ret_val + poly->coeffs[i];
    double xmul = 1.0;
    for(int j = 0; j < (poly->power[i] - poly->power[i-1]); j++){
      xmul = xmul * x;
    }
    ret_val = ret_val * xmul;	  
  }
  ret_val = ret_val + poly->coeffs[0];
  
  for(int j = 0; j < poly->power[0]; j++){
    ret_val = ret_val * x;
  }  
  return ret_val;
#endif  

}

