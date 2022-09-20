#include <math.h>
#include "pi.h"
#include "rlibm.h"
#include "sinpicospi.h"

#define PI_OVER_512 0x1.921fb54442d18p-8 

double rlibm_sinf(float x) {
  float_x fx;
  fx.f = x; 
  int s = (fx.x & 0x80000000) ? -1 : 1;
  fx.x &= 0x7FFFFFFF;

  if (fx.x >= 0x7F800000) {
    return 0.0f/0.0f;
  }
  if (fx.x == 0x4AFDECE4) {
    return -0x1.ff6dc18p-1*s;
  }
  if (fx.x == 0x55CAFB2A) {
    return -0x1.fcf42d8p-1*s;
  }
  if (fx.x == 0x73243F06) {
    return 0x1.2875078p-2*s;
  }
  double y;
  if (fx.x < 0x39E89769) {
    y = fx.f*0x1.ffffffffffffep-1; 
  } else {
    double R;
    int64_t N;

    if (fx.x <= 0x48000000) {
      double prod = fx.f * _512_over_pi_28[0];
      double kd = floor(prod);
      N = (int64_t) kd;
      R = prod - kd;
      R = fma(fx.f, _512_over_pi_28[1], R);
      R = fma(fx.f, _512_over_pi_28[2], R);
      if (R > 1.0f) {
	N += 1;
	R -= 1.0f;
      }
    } else {
      int x_biased_e = fx.x >> 23L;
      int x_lsb_e = x_biased_e - 150;
      int idx = 0;
      while (x_lsb_e + _512_over_pi_28_exp[idx] > 0) {
	idx++;
      }
      double prod_hi =  (idx) ? fx.f * _512_over_pi_28[idx-1] : 0;
      int64_t k_hi_int = (int64_t) prod_hi;
      k_hi_int &= 0x3FF;
      double k_hi = (double) k_hi_int;
      double prod_mid = fx.f * _512_over_pi_28[idx];
      double k_mid = floor(prod_mid);
      double frac = prod_mid - k_mid;
      double prod_lo = fma(fx.f, _512_over_pi_28[idx + 1], frac);
      double k_lo = floor(prod_lo); 
      N = (int64_t)(k_hi + k_mid + k_lo);

      R = prod_lo - k_lo;
      R = fma(fx.f, _512_over_pi_28[idx + 2], R);
      R = fma(fx.f, _512_over_pi_28[idx + 3], R); 
      if (R > 1.0f) {
	N += 1;
	R -= 1.0f;
      }
    }
    R *= PI_OVER_512;
    
    unsigned N2 = N & 0xFF;
    unsigned I = N >> 8;
    
    if (I & 0x1) {
      N2 = 255 - N2;
      R = PI_OVER_512 - R;
    } 
   
    if (I & 0x2) s *= -1;
   
    double R2 = R * R;

    double sinpiM = sinpiMBy512[N2];
    double cospiM = cospiMBy512[N2];
    double sinpiR, cospiR;

    sinpiR = 0x1.110dd89739b5cp-7;
    sinpiR *= R2;
    sinpiR += -0x1.55555554ce812p-3;
    sinpiR *= R2;
    sinpiR += 0x1.ffffffffffffep-1;
    sinpiR *= R;

    cospiR = 0x1.5553dd3610b9p-5;
    cospiR *= R2;
    cospiR += -0x1.ffffffff9c28fp-2; 
    cospiR *= R2;
    cospiR += 0x1.fffffffffffffp-1; 

    y = sinpiM * cospiR + cospiM * sinpiR;
  }
  y *= s;
  return y;
}
