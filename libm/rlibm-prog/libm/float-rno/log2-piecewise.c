#include "float-rno/float_rno_lib.h"
#include "float_lib.h"
#include "float-rno/log2.h"
#include <math.h>

double rlibm_rno_log2_piecewise(float x) {
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

  double y;
  
  if (f < 3.906249999999999566319131005798226397018879652023315429687500e-03) {
    y = -3.5668687621348926786168931357678957283496856689453125000000000000000000e-01;
    y *= f;
    y += 4.8087976044912084105931171507108956575393676757812500000000000000000000e-01;
    y *= f;
    y += -7.2134748569212248092696881940355524420738220214843750000000000000000000e-01;
    y *= f;
    y += 1.4426950408671799230830856686225160956382751464843750000000000000000000e+00;
    y *= f;
  } else {
    y = -3.5339430743329897088855773290561046451330184936523437500000000000000000e-01;
    y *= f;
    y += 4.8083087641098110065485116138006560504436492919921875000000000000000000e-01;
    y *= f;
    y += -7.2134725011798661586936987077933736145496368408203125000000000000000000e-01;
    y *= f;
    y += 1.4426950404954139717261796249658800661563873291015625000000000000000000e+00;
    y *= f;
  }
    
  y += __log2_lut[FIndex];
  y += m;
  
  return y;
}
