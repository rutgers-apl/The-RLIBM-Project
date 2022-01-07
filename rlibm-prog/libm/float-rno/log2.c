#include "float-rno/float_rno_lib.h"
#include "float_lib.h"
#include "float-rno/log2.h"
#include <math.h>

double rlibm_rno_log2(float x) {
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

#define C1 1.4426950408932217584379031904973089694976806640625000000000000000000000e+00
#define C2 -7.2134752833423521067857109301257878541946411132812500000000000000000000e-01
#define C3 4.8090209033186281928351490932982414960861206054687500000000000000000000e-01
#define C4 -3.6134059092802622847884208567847963422536849975585937500000000000000000e-01 
#define C5 3.2718656156289893655042533282539807260036468505859375000000000000000000e-01

  double y = 0.0;

  y = C5;
  y *= f;
  y += C4;
  y *= f;
  y += C3;
  y *= f;
  y += C2;
  y *= f;
  y += C1;
  y *= f;
    
  y += __log2_lut[FIndex];
  y += m;
  
  return y;
}
