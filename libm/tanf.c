/*
MIT License

Copyright (c) 2022 Sehyeok Park and Santosh Nagarakatte, Rutgers
Architecture and Programming Languages (RAPL) Group

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/


#include <math.h>
#include "pi.h"
#include "rlibm.h"
#include "sinpicospi.h"

#define PI_OVER_512 0x1.921fb54442d18p-8 

double rlibm_tanf(float x) {
  float_x fx;
  fx.f = x; 
  int s = (fx.x & 0x80000000) ? -1 : 1;
  fx.x &= 0x7FFFFFFF;
  if (fx.x >= 0x7F800000) {
    return 0.0f/0.0f;
  }
  if (fx.x == 0x3F8A1F62) {
    return 0x1.ddf9f58p+0*s;
  }
  if (fx.x == 0x4D56D355) {
    return 0x1.e803048p-3*s;
  }
  if (fx.x == 0x5FFD33A4) {
    return 0x1.a0d9178p+0*s;
  }
  if (fx.x == 0x6AD36709) {
    return -0x1.c5612e8p-1*s;
  }
  double y;
  if (fx.x < 0x39B89BA3) {
    y = fx.f * 0x1.0000008000007p+0*s;
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
   
    double sinpiM, cospiM, sinpiR, cospiR;
    unsigned N2, I;

    double sinR = R;

    N2 = N & 0xFF;
    I = N >> 8;
  	  
    if (I & 0x1) {
      N2 = 255 - N2;
      sinR = PI_OVER_512 - sinR;
    } 
   
    if (I & 0x2) s *= -1;
   
    double sinR2 = sinR * sinR;

    sinpiM = sinpiMBy512[N2];
    cospiM = cospiMBy512[N2];
    sinpiR, cospiR;

    sinpiR = 0x1.110dd89739b5cp-7;
    sinpiR *= sinR2;
    sinpiR += -0x1.55555554ce812p-3;
    sinpiR *= sinR2;
    sinpiR += 0x1.ffffffffffffep-1;
    sinpiR *= sinR;

    cospiR = 0x1.5553dd3610b9p-5;
    cospiR *= sinR2;
    cospiR += -0x1.ffffffff9c28fp-2; 
    cospiR *= sinR2;
    cospiR += 0x1.fffffffffffffp-1; 
    
    double sin_x = sinpiM * cospiR + cospiM * sinpiR;
    sin_x *= s;

    double cosR = R;

    N2 = N & 0xFF;
    I = (N >> 8) + 1;
    s = (I & 0x2) ? -1 : 1;

    if (I & 1) {
      if (N2 == 0) {
	cospiM = 1.0;
	sinpiM = 0.0;
      } else {
	N2++;
	cosR = PI_OVER_512 - cosR;
	cospiM = sinpiMBy512[256 - N2];
	sinpiM = cospiMBy512[256 - N2];
      }
    } else {
      cospiM = sinpiMBy512[N2];
      sinpiM = cospiMBy512[N2];
    }

    double cosR2 = cosR * cosR;

    sinpiR = 0x1.110dd89739b5cp-7;
    sinpiR *= cosR2;
    sinpiR += -0x1.55555554ce812p-3;
    sinpiR *= cosR2;
    sinpiR += 0x1.ffffffffffffep-1;
    sinpiR *= cosR;

    cospiR = 0x1.5553dd3610b9p-5;
    cospiR *= cosR2;
    cospiR += -0x1.ffffffff9c28fp-2; 
    cospiR *= cosR2;
    cospiR += 0x1.fffffffffffffp-1; 

    double cos_x = sinpiM * sinpiR + cospiM * cospiR;
    cos_x *= s;

    y = sin_x/cos_x;
  }
  return y;
}
