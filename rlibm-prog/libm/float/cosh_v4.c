#include "float_lib.h"
#include "sinhcosh.h"

#define CONST64BYLN2 92.332482616893656768297660164535045623779296875
#define LN2BY64 0.01083042469624914509729318723429969395510852336883544921875

float rlibm_fast_cosh_v4(float x) {
  float_x fx;
  fx.f = x;
  fx.x &= 0x7FFFFFFF;
  
  // Take care of special cases
  if (fx.x <= 968164595) return 1.0;
  
  if (fx.x >= 1119016189) {
      if (fx.x > 0x7F800000) return 0.0f/0.0f;
      return 1.0f / 0.0f;
  }
  
  // Perform range reduction
  double xp = fx.f * CONST64BYLN2;
  int N = (int)xp;
  int N2 = N % 64;
  if (N2 < 0) N2 += 64;
  int N1 = N - N2;
  int I = N1 / 64;
  // redued input x'
  double R = fx.f - N * LN2BY64;
  // x' squared
  double R2 = R * R;
  
  double sinhHigh = sinhKLn2[I];
  double coshHigh = coshKLn2[I];
  double sinhMid = sinhKLn2By64[N2];
  double coshMid = coshKLn2By64[N2];
  
  double sinhHM = sinhHigh * coshMid + coshHigh * sinhMid;
  double coshHM = sinhHigh * sinhMid + coshHigh * coshMid;
  
  // Compute sinh component
  double sinhL;
  sinhL = 4.919193754695290465850376904199947603046894073486328125e-02;
  sinhL *= R2;
  sinhL += 8.32675387692339251388684573385035037063062191009521484375e-03;
  sinhL *= R2;
  sinhL += 1.666666668375217097430862622786662541329860687255859375e-01;
  sinhL *= R2;
  sinhL += 9.999999999999997779553950749686919152736663818359375e-01;
  sinhL *= R;
  
  // Compute cosh component
  double coshL;
  coshL = 4.919193754695290465850376904199947603046894073486328125e-02;
  coshL *= R2;
  coshL += 8.32675387692339251388684573385035037063062191009521484375e-03;
  coshL *= R2;
  coshL += 1.666666668375217097430862622786662541329860687255859375e-01;
  coshL *= R2;
  coshL += 9.999999999999997779553950749686919152736663818359375e-01;
  
  
  // Perform output compensation
  return sinhHM * sinhL + coshHM * coshL;
}
