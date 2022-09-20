#include "float-rno/float_rno_lib.h"
#include "float_lib.h"
#include "float-rno/sinpicospi.h"

#define PI 3.141592653589793115997963468544185161590576171875

double rlibm_rno_sinpi(float x) {
  float_x fX;
  fX.f = x;
  unsigned long s = (fX.x & 0x80000000) == 0 ? 0x0 : 0x8000000000000000;
  fX.x &= 0x7FFFFFFF;
  
  if (fX.x <= 872158519) {
    if (fX.x == 869073367) {
      if (s) return -2.9990871297513876925222575664520263671875e-07;
      return 2.9990871297513876925222575664520263671875e-07;
    }

    return (double)x * PI;
  }
  
  if (fX.x >= 0x4b000000) {
    if (fX.x >= 0x7F800000) {
        return 0.0/0.0;
    }
    return 0.0;
  }
  
  double xp = fX.f * 512;
  unsigned N = (unsigned)xp;
  unsigned N2 = N & 0xFF;
  unsigned I = N >> 8;
  double R = fX.f - N * 0.001953125;
  if (I & 0x2) s ^= 0x8000000000000000;
  
  double_x dX;
    
  // More special cases when sinpi(x) is exactly representable
  if (R == 0 && N2 == 0) {
    if (I & 0x1) {
      dX.d = 1.0;
      dX.x |= s;
      return dX.d;
    } else {
      return 0.0;
    }
  }

  if (I & 0x1) {
      N2 = 255 - N2;
      R = 0.001953125 - R;
  }

  dX.d = R;
  double R2 = R * R;
  
  double cospiR, sinpiR;
  
  sinpiR = 2.552988167646617068129444305668585002422332763671875;
  sinpiR *= R2;
  sinpiR += -5.16771278732025063362698347191326320171356201171875;
  sinpiR *= R2;
  sinpiR += 3.141592653589793115997963468544185161590576171875;
  sinpiR *= R;
  
  cospiR = 4.060478460719021853719823411665856838226318359375;
  cospiR *= R2;
  cospiR += -4.934802206516291533944240654818713665008544921875;
  cospiR *= R2;
  cospiR += 1.0000000000000011102230246251565404236316680908203125;
      
  dX.d = sinpiMBy512[N2] * cospiR + cospiMBy512[N2] * sinpiR;
  dX.x ^= s;
  
  return dX.d;
}
