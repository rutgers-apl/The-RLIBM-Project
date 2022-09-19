#include <math.h>
#include "rlibm.h"
#include "exp2.h"

double rlibm_exp2(float x) {
  float_x fx;
  fx.f = x;
  
  // Take care of special cases
  if ((fx.x & 0x7FFFFFFF) == 0) return 1.0;
  
  if (fx.x <= 876128826) {
    if (fx.x <= 867740218) return 1.0000000298023223876953125;
    return 1.0000000894069671630859375;
  }
  
  if (1124073472 <= fx.x && fx.x <= 3015223867) {
    if (fx.x < 0x80000000) {
      if (fx.x < 0x7F800000) return 3.40282361850336062550457001444955389952e+38;
      if (fx.x == 0x7F800000) return 1.0 / 0.0;
      return 0.0/0.0;
    }
    
    // negative counterpart
    if (fx.x <= 3006835259) return 0.99999998509883880615234375;
    
    return 0.99999995529651641845703125;
  }
  
  if (fx.x >= 3272998913) {
    if (fx.x == 0xFF800000) return 0.0;
    if (fx.x < 0xFF800000) return ldexp(1.0, -151);
    return 0.0/0.0;
  }
  
  // Perform range reduction
  double xp = x * 64;
  int N = (int)xp;
  int N2 = N % 64;
  if (N2 < 0) N2 += 64;
  int N1 = N - N2;
  
  int M = N1 / 64;
  int J = N2;
  double R = x - N * 0.015625;
  
  if (R == 0.0 && N2 == 0) {
    return ldexp(1.0, M);
  }

  double y = 0.0;

  if(R == 0x1.853a6ep-9){
    y = 0x1.0087090000001p+0 ;
  }
  else{
    double coeffs[] = {
		       0x1.fffffffffffffp-1,
		       0x1.62e42fefa3806p-1,
		       0x1.ebfbdff81263ap-3,
		       0x1.c6b08b875921bp-5,
		       0x1.3b2b5ff31806dp-7,
		       0x1.62bcac112f1e1p-10,
    };

    double xsquare = R * R;
    double temp1 = fma(R, coeffs[1], coeffs[0]);
    double temp2 = fma(R, coeffs[5], coeffs[4]);
    double temp3 = fma(R, coeffs[3], coeffs[2]);

    double temp4 = fma(xsquare, temp2, temp3);
    y = fma(xsquare, temp4, temp1);
  }
  
  // Perform output compensation
  return y * ldexp(exp2JBy64[J], M);
}
