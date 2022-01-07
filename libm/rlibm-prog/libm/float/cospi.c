#include "float_lib.h"
#include "sinpicospi.h"

float rlibm_fast_cospi(float x) {
  float_x fX;
  fX.f = x;
  fX.x &= 0x7FFFFFFF;
  
  // Special cases:
  // If x is smaller than 0x38a2f983, then it's 1.0f
  if (fX.x <= 0x38a2f983) {
      return 1.0f;
  }
  
  if (fX.x >= 0x4b000000) {
      // If x >= 0x7F800000, then result is NaN
      if (fX.x >= 0x7F800000) return 0.0f/0.0f;
      // If x >= 2^24, then result is always 1.0f
      if (fX.x >= 0x4b800000) return 1.0f;
      // If x >= 2^23, then if x is even, then 1.0f
      if ((fX.x & 0x1) == 0) return 1.0f;
      // Otherwise, then -1.0f
      return -1.0f;
  }
  
  // Range Reduction
  double xp = fX.f * 512.0;
  unsigned N = (unsigned)xp;
  unsigned N2 = N & 0xFF;
  unsigned I = (N >> 8) + 1;
  double R, cospiM, sinpiM;
  unsigned s = (I & 0x2) ? 0x80000000 : 0;

  if (I & 1) {
    if (N2 == 0) {
      R = fX.f - N * 0.001953125;
      cospiM = 1.0;
      sinpiM = 0.0;
    }
    else {
      N2++;
      R = (N + 1) * 0.001953125 - fX.f;
      cospiM = sinpiMBy512[256 - N2];
      sinpiM = cospiMBy512[256 - N2];
    }
  } else {
    R = fX.f - N * 0.001953125;
    cospiM = sinpiMBy512[N2];
    sinpiM = cospiMBy512[N2];
  }
  
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

  fX.f = cospiM * cospiR + sinpiM * sinpiR;
  fX.x ^= s;
  
  return fX.f;
}
