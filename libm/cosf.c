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

double rlibm_cosf(float x) {
  float_x fx;
  fx.f = x; 
  fx.x &= 0x7FFFFFFF;
  if (fx.x >= 0x7F800000) {
    return 0.0f/0.0f;
  }
  if (fx.x == 0) {
    return 1.0;
  }
  if (fx.x == 0x491A2430) {
    return -0x1.edfe2f8p-1;
  }
  if (fx.x == 0x55325019) {
    return 0x1.9d4ba48p-1;
  }
  if (fx.x == 0x59443C0A) {
    return 0x1.84bec48p-1;
  }
  if (fx.x == 0x5AA4542C) {
    return 0x1.f481488p-2;
  }
  if (fx.x == 0x6115CB11) {
    return 0x1.f0285d8p-1;
  }
  if (fx.x == 0x61703976) {
    return 0x1.b598ab8p-2;
  }
  if (fx.x == 0x6F5D2D4C) {
    return -0x1.ac093c8p-1;
  }

  int s = 1;
  double y;
  if (fx.x < 0x39800001) {
    y = 0x1.fffffffffffffp-1;
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

    double sinpiM;
    double cospiM;
    
    unsigned N2 = N & 0xFF;
    unsigned I = (N >> 8) + 1;
    s = (I & 0x2) ? -1 : 1;

    if (I & 1) {
      if (N2 == 0) {
	cospiM = 1.0;
	sinpiM = 0.0;
      } else {
	N2++;
	R = PI_OVER_512 - R;
	cospiM = sinpiMBy512[256 - N2];
	sinpiM = cospiMBy512[256 - N2];
      }
    } else {
      cospiM = sinpiMBy512[N2];
      sinpiM = cospiMBy512[N2];
    }
   
    double R2 = R * R;

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
    
    y = sinpiM * sinpiR + cospiM * cospiR;
  }
  y *= s;
  return y;
}
