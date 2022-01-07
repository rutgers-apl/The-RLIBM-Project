#include "float-rno/float_rno_lib.h"
#include "float_lib.h"
#include "float-rno/sinpicospi.h"

double rlibm_rno_cospi(float x) {
  float_x fX;
  fX.f = x;
  fX.x &= 0x7FFFFFFF;

  if (fX.x <= 950204803) {
    if (fX.x == 0) return 1.0;
    return 0.99999998509883880615234375;
  }
  
  if (fX.x >= 0x4b000000) {
    if (fX.x >= 0x7F800000) return 0.0 / 0.0;
    if (fX.x >= 0x4b800000) return 1.0;
    
    // If x >= 2^23, then if x is even, then 1.0f
    if ((fX.x & 0x1) == 0) return 1.0;
    
    // Otherwise, then -1.0f
    return -1.0;
  }
    
  // Range Reduction
  double xp = fX.f * 512.0;
  unsigned N = (unsigned)xp;
  unsigned N2 = N & 0xFF;
  unsigned I = (N >> 8) + 1;
  double R, cospiM, sinpiM;
  unsigned long s = (I & 0x2) ? 0x8000000000000000 : 0;
  R = fX.f - N * 0.001953125;
  
  if (R == 0 && N2 == 0) {
    if (I & 1l) {
      double_x resx;
      resx.d = 1.0;
      resx.x |= s;
      return resx.d;
    } else {
      return 0.0;
    }
  }
    
  if (I & 1) {
    if (N2 == 0) {
      cospiM = 1.0;
      sinpiM = 0.0;
    }
    else {
      N2++;
      R = 0.001953125 - R;
      cospiM = sinpiMBy512[256 - N2];
      sinpiM = cospiMBy512[256 - N2];
    }
  } else {
    cospiM = sinpiMBy512[N2];
    sinpiM = cospiMBy512[N2];
  }

  double R2 = R * R;
  double cospiR, sinpiR;
  
  sinpiR = 2.5550488062509071340855371090583503246307373046875;
  sinpiR *= R2;
  sinpiR += -5.16771279405369199366759858094155788421630859375;
  sinpiR *= R2;
  sinpiR += 3.14159265358979755689006196917034685611724853515625;
  sinpiR *= R;
  
  cospiR = 4.0601744557652210687592742033302783966064453125;
  cospiR *= R2;
  cospiR += -4.93480220469346253509002053760923445224761962890625;
  cospiR *= R2;
  cospiR += 1.0000000000000015543122344752191565930843353271484375;

  double_x dX;
  dX.d = cospiM * cospiR + sinpiM * sinpiR;
  dX.x ^= s;
  
  return dX.d;
}
