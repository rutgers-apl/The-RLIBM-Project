#include "polygen.h"

double rlibm_poly_evaluation(double x, polynomial* poly){

  if(poly->termsize == 1){
      return x*poly->coeffs[0];
  } else if(poly->termsize == 2){
    double temp = x * x * poly->coeffs[1];
    return fma(x, poly->coeffs[0],temp);
    
    //      return x*poly->coeffs[0] + x*x*poly->coeffs[1];
  } else if(poly->termsize == 3){
    double temp = x * x * fma(x, poly->coeffs[2], poly->coeffs[1]);
    return fma(x, poly->coeffs[0], temp);
    
    //     return x*poly->coeffs[0] + x*x*(poly->coeffs[1] + x*poly->coeffs[2]);
  } else if(poly->termsize == 4){

    double temp1 = fma(x, poly->coeffs[2], poly->coeffs[1]);
    double xsquare = x*x;
    double temp2 = fma(xsquare, poly->coeffs[3], temp1);
    double temp3 = xsquare * temp2;
    return fma(x, poly->coeffs[0], temp3);
    
    //      return x*poly->coeffs[0] + x*x*(poly->coeffs[1] + x*poly->coeffs[2] + x*x*poly->coeffs[3]);
  } else if(poly->termsize == 5){
    double xsquare = x*x;
    double temp1 = fma(x, poly->coeffs[4], poly->coeffs[3]);

    double temp2 = fma(x, poly->coeffs[2], poly->coeffs[1]);
    double temp3 = fma(xsquare, temp1, temp2);
    double temp4 = xsquare * temp3;
    return fma(x, poly->coeffs[0], temp4);
    
    //     return x*poly->coeffs[0] + x*x*(poly->coeffs[1] + x*poly->coeffs[2] + x*x*(poly->coeffs[3] + x*poly->coeffs[4]));

    
  } else if(poly->termsize == 6){
    double xsquare = x* x;
    double temp1 = fma(x, poly->coeffs[4], poly->coeffs[3]);
    double temp2 = fma(xsquare, poly->coeffs[5], temp1);

    double temp3 = fma(x, poly->coeffs[2], poly->coeffs[1]);
    double temp4 = fma(xsquare, temp2, temp3);
    double temp5 = xsquare * temp4;
    return fma(x, poly->coeffs[0], temp5);

    
    //   return x*poly->coeffs[0] + x*x*(poly->coeffs[1] + x*poly->coeffs[2] + x*x*(poly->coeffs[3] + x*poly->coeffs[4] + x*x*poly->coeffs[5]));
  } else if(poly->termsize == 7){
    double xsquare = x * x;
    double temp1 = fma(x, poly->coeffs[6], poly->coeffs[5]);
    double temp2 = fma(x, poly->coeffs[4], poly->coeffs[3]);
    double temp3 = fma(xsquare, temp1, temp2);
    double temp4 = fma(x, poly->coeffs[2], poly->coeffs[1]);

    double temp5 = xsquare * (temp4 + temp3);
    return fma(x, poly->coeffs[0], temp5);
    
    //    return x*poly->coeffs[0] + x*x*(poly->coeffs[1] + x*poly->coeffs[2] + x*x*(poly->coeffs[3] + x*poly->coeffs[4] + x*x*(poly->coeffs[5] + x*poly->coeffs[6])));
  }

  // default horner's method
  
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
}

