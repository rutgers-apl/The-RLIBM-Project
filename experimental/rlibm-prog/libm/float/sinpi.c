#include "float_lib.h"
#include "sinpicospi.h"

#define PI 3.141592653589793115997963468544185161590576171875

float rlibm_fast_sinpi(float x) {
  float_x fX;
  fX.f = x;
  unsigned s = fX.x & 0x80000000;
  fX.x &= 0x7FFFFFFF;
  
  // Special cases:
  if (fX.x <= 0x33fc1537) {
    return PI * (double)x;
  }
  
  if (fX.x >= 0x4b000000) {
    if (fX.x >= 0x7F800000) {
      return 0.0f/0.0f;
    }
    return 0.0f;
  }
  
  double xp = fX.f * 512;
  unsigned N = (unsigned)xp;
  unsigned N2 = N & 0xFF;
  unsigned I = N >> 8;
  double R;
  
  if (I & 0x1) {
    N2 = 255 - N2;
    R = (N + 1) * 0.001953125 - fX.f;
  } else R = fX.f - N * 0.001953125;
  
  if (I & 0x2) s ^= 0x80000000;
  
  double R2 = R * R;
  double cospiR, sinpiR;
  sinpiR = 2.553984001513528223625826285569928586483001708984375;
  sinpiR *= R2;
  sinpiR += -5.16771279048179632553683404694311320781707763671875;
  sinpiR *= R2;
  sinpiR += 3.141592653589793560087173318606801331043243408203125;
  sinpiR *= R;
  
  cospiR = 4.0605714171214888352778871194459497928619384765625;
  cospiR *= R2;
  cospiR += -4.93480220691509430253063328564167022705078125;
  cospiR *= R2;
  cospiR += 1.0000000000000028865798640254070051014423370361328125;
        
  fX.f = sinpiMBy512[N2] * cospiR + cospiMBy512[N2] * sinpiR;
  fX.x ^= s;
  
  return fX.f;
}
