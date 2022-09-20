#include "rlibm.h"
#include "log2.h"
#include <math.h>

double rlibm_log2(float x) {
  float_x fix, fit, spec;
  fix.f = x;
  int m = 0;

  if (fix.x < 0x800000 || fix.x >= 0x7F800000) {
    if ((fix.x & 0x7FFFFFFF) == 0) { // log(+/-0) = -infty
      fix.x = 0xFF800000;
      return fix.f;
    }
    
    if (fix.x > 0x7FFFFFFF) { // Log(-val) = NaN
      return (x - x) / 0;
        
    }
    
    if (fix.x >= 0x7F800000) {
      return x + x;
    }

    // Special case when we have denormal input and exact result
    int exp;
    spec.f = frexpf(fix.f, &exp);
    if (spec.x == 0x3f000000) return (double)(exp - 1);

    fix.f *= 8.388608e+06;
    m -= 23;
  }
  
  m += fix.x >> 23;
  m -= 127;
  fix.x &= 0x007FFFFF;
  
  if (fix.x == 0) {
    return (double)m;
  }
  
  fix.x |= 0x3F800000;
  
  fit.x = fix.x & 0x007F0000;
  int FIndex = fit.x >> 16;
  fit.x |= 0x3F800000;
  
  double f = fix.f - fit.f;
  f *= __log_oneByF[FIndex];

  double y = 0.0;


  // FMA coefficients
  double coeffs[] = {
		     1.4426950408932206482148785653407685458660125732421875000000000000000000e+00,
		     -7.2134752833238291458428648184053599834442138671875000000000000000000000e-01,
		     4.8090208952661678276641055163054261356592178344726562500000000000000000e-01,
		     -3.6134046214924120388189976438297890126705169677734375000000000000000000e-01,
		     3.2717964062254517587646773790766019374132156372070312500000000000000000e-01		     
  };
  
  double xsquare = f*f;
  double temp1 = fma(f, coeffs[4], coeffs[3]);
  double temp2 = fma(f, coeffs[2], coeffs[1]);
  double temp3 = fma(xsquare, temp1, temp2);
  double temp4 = xsquare * temp3;
  y = fma(f, coeffs[0], temp4);    
  y += __log2_lut[FIndex];
  y += m;
  
  return y;
}
