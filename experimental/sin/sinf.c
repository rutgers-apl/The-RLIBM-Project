
/* MIT License */

/* Copyright (c) 2022 Sehyeok Park and Santosh Nagarakatte, The RLIBM
   project, Rutgers Architecture and Programming Languages (RAPL)
   Group */

/* Permission is hereby granted, free of charge, to any person obtaining a copy */
/* of this software and associated documentation files (the "Software"), to deal */
/* in the Software without restriction, including without limitation the rights */
/* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
/* copies of the Software, and to permit persons to whom the Software is */
/* furnished to do so, subject to the following conditions: */

/* The above copyright notice and this permission notice shall be included in all */
/* copies or substantial portions of the Software. */

/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
/* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
/* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
/* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
/* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
/* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
/* SOFTWARE. 
*/


#include "sinf.h"


float rlibm_fast_sin(float x) {
  FLT fX;
  fX.f = x;

  /* look at the sign */
  int s = (fX.i & 0x80000000) ? -1 : 1;
  /* fX is the absolute value */
  fX.i &= 0x7FFFFFFF;

  /* handle NANs */
  if (fX.i >= 0x7F800000) {
    return 0.0f/0.0f;
  }

  /* for small values of x up to  0x1.d12ed2p-12, sin(x) = x */
  if (fX.i < 0x39E89769) {
    return fX.f*s;
  }

  /* one special case input: +/- 0x1.487e0cp+103  */
  if (fX.i == 0x73243F06) {
    fX.i = 0x3E943A84;
    return fX.f*s;
  }

  /* Range reduction happens at two levels:
     First level:
     express input x = k*pi + r, k is the nearest integer.
     sin(x) = sin(k*pi + r) = sin(k*pi)cos(r) + cos(k*pi)sin(r). 
     sin(x) = (-1)^k. sin(r).

     To compute sin(r), we use sinpi(y), where y is x/pi -k.
     r = x - k*pi = pi( x/pi - k).
     sin(r) = sinpi(y).
     Here |y| < 0.5

     Second level:  
     For sinpi(y), we use the range reduction in Section 2 of our PLDI 2021 paper:
     https://people.cs.rutgers.edu/~sn349/papers/rlibm32-pldi-2021-preprint.pdf

     sinpi(x) = (-1)^K1 (sinpi(N/512) cospi(R) + cospi(N/512) sinpi(R))
     Here K1 and R are computed as described in our above paper.
  */
  double R;
  /* Get r */
  if (fX.i <= 0x4C800000) {
    /* small range reduction */
    double prod = fX.f * one_over_pi_28[0];
    double kd = round(prod);
    /* max value of kd is less than the range of an int */
    int n = ((int64_t)kd) & 0x1;

    R = prod - kd;
    R = fma(fX.f, one_over_pi_28[1], R);
    R = fma(fX.f, one_over_pi_28[2], R);
    n ^= (R < 0);
    if (R > 0.5) {
      R = 1.0f - R;
    }
    R = fabs(R);
    if (n & 0x1) s *= -1;
  } else {
    
    int x_biased_e = (fX.i & 0x7F800000) >> 23L;
    int x_lsb_e = x_biased_e - 150;
    int idx = 0;

    /* we are only interested in the last bit of k, ignore all other
       bits */
  
    while (x_lsb_e + one_over_pi_28_exp[idx] > 0) {
      idx++;
    }
    double prod_hi = fX.f * one_over_pi_28[idx];

    /* max value of k_hi is less than the range of an int */
    double k_hi = round(prod_hi);
    double frac = prod_hi - k_hi;
    double prod_lo = fma(fX.f, one_over_pi_28[idx + 1], frac);
    double k_lo = round(prod_lo);
    int64_t k = (int64_t)(k_hi + k_lo);
    int n = k & 0x1;

    R = prod_lo - k_lo;
    R = fma(fX.f, one_over_pi_28[idx + 2], R);
    R = fma(fX.f, one_over_pi_28[idx + 3], R); 
    n ^= (R < 0);
    if (R > 0.5) {
      R = 1.0f - R;
    }
    R = fabs(R);
    if (n & 0x1) s *= -1;
  }

  double R512 = R * 512;
  unsigned N = (unsigned) R512;
  unsigned N2 = N & 0xFF;
  R -= N*0.001953125;

  double R2 = R * R;
  double cospiR, sinpiR;

  sinpiR = 0x1.466a82b601746p+1;
  sinpiR *= R2;
  sinpiR += -0x1.4abbce626862bp+2;
  sinpiR *= R2;
  sinpiR += 0x1.921fb54442d19p+1;
  sinpiR *= R;

  cospiR = 0x1.03c0a464c9a04p+2;
  cospiR *= R2;
  cospiR += -0x1.3bd3cc9bb1ac6p+2; 
  cospiR *= R2;
  cospiR += 0x1p+0;
  
  double y = sinpiMBy512[N2] * cospiR + cospiMBy512[N2] * sinpiR;
  fX.f = y;
  fX.f *= s;
  return fX.f;
}

