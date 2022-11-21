/*
MIT License

Copyright (c) 2022 Sehyeok Park, Mridul Aanjaneya, and Santosh
Nagarakatte, Rutgers Architecture and Programming Languages (RAPL)
Group

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
#include "rlibm.h"
#include "sinpicospi.h"
#include "pi.h"


#define PI_OVER_256 0x1.921fb54442d18p-7 

typedef unsigned __int128 uint128_t;


double rlibm_sinf(float x) {
  double_x dY;
  float_x fX;
  fX.f = x;
  fX.x &= 0x7FFFFFFF;
  uint64_t sgn = (x < 0) ? 0x8000000000000000 : 0;
  double R, R2;
  if (fX.x < 0x39E89769) {
    dY.d = fX.f*0x1.ffffffffffffep-1;
    dY.x ^= sgn;
    return dY.d;
  } else if (fX.x < 0x7F800000) {	  
    int s = (fX.x >> 23) - 150;
    uint64_t m = (fX.x&0x7FFFFF)|1<<23;
    int N;
    uint64_t a;
    if (fX.x < 0x4A000000) {
      if (fX.x < 0x3CC90FDB) {
	R = fX.f;
	R2 = R*R;
	double temp = fma(R2, 0x1.11130658ac4e4p-7, -0x1.55555558c07f5p-3);
	dY.d = R*fma(temp, R2, 0x1.0000000000011p+0);
	dY.x ^= sgn;
	return dY.d;
      }
      if (fX.x == 0x3F3ADC51) {
	dY.d = 0x1.5568898p-1;
	dY.x ^= sgn;
	return dY.d;
      }
      if (fX.x == 0x4371ADE3) {
	dY.d = 0x1.c5b4ab8p-3;
	dY.x ^= sgn;
	return dY.d;
      }

      /* 80 bits of 1/pi stored as 40-bit pieces in two 64-bit
	 integers */
      uint64_t p0 = m*0x441529fc27;
      uint64_t p1 = m*0xa2f9836e4e; p1+=(p0>>40);
      N = (p1>>(33-s));
      a = p1<<(31+s)|((p0<<24)>>(33-s));
    } else {
      if (fX.x == 0x55CAFB2A) {
	dY.d = -0x1.fcf42d8p-1;
	dY.x ^= sgn;
	return dY.d;
      }
      if (fX.x == 0x5DADD689) {
	dY.d = -0x1.e9f93b8p-1;
	dY.x ^= sgn;
	return dY.d;
      }
      if (fX.x == 0x67A9242B) {
	dY.d = -0x1.ff57018p-1;
	dY.x ^= sgn;
	return dY.d;
      }
      if (fX.x == 0x73243F06) {
	dY.d = 0x1.2875078p-2;
	dY.x ^= sgn;
	return dY.d;
      }

      /* 256 bits of 1/pi for range reduction with large arguments */
      static const uint64_t ipi[] = {0xfe5163abdebbc562, 0xdb6295993c439041, 0xfc2757d1f534ddc0, 0xa2f9836e4e441529};
      uint128_t p0 = (uint128_t)m*ipi[0];
      uint128_t p1 = (uint128_t)m*ipi[1]; p1 += p0>>64;
      uint128_t p2 = (uint128_t)m*ipi[2]; p2 += p1>>64;
      uint128_t p3 = (uint128_t)m*ipi[3]; p3 += p2>>64;
      uint64_t p3h = p3>>64, p3l = p3, p2l = p2, p1l = p1;
      if (s < 57) {
	N = (p3h<<(7+s))|(p3l>>(57-s));
	a = (p3l<<(7+s))|(p2l>>(57-s));
      } else if (s == 57) {
	N = p3l;
	a = p2l;
      } else {
	N = (p3l<<(s-57))|(p2l>>(121-s));
	a = (p2l<<(s-57))|(p1l>>(121-s));
      }
    } 
    long N2 = N & 0x7F;
    unsigned I = N >> 7;
    if (I & 0x2) sgn ^= 0x8000000000000000;
    N2 ^= ~(0xFFFFFFFF+(I & 0x1));
    N2 += 128&((I & 0x1)<<7);
    a ^= ~(0xFFFFFFFFFFFFFFFF+(I & 0x1));
    R = a*0x1p-64*PI_OVER_256;

    double sinpiM = sinpiMBy256[(unsigned)N2];
    double cospiM = cospiMBy256[(unsigned)N2];
 
    double temp1, temp2;
    R2 = R*R;

    temp1 = fma(R2, 0x1.11154bcb5c9e2p-7, -0x1.555555583bf78p-3);
    dY.d = cospiM*R*fma(temp1, R2, 0x1.000000000000cp+0);
    if (N2 == 0) {
      dY.x ^= sgn;
      return dY.d;
    }
    temp2 = fma(R2, 0x1.5554a56a7705fp-5, -0x1.ffffffffad9f2p-2);
    dY.d += sinpiM*fma(temp2, R2, 0x1.ffffffffffffep-1);
    
    dY.x ^= sgn;
    return dY.d;
  } else { 
    return 0.0f/0.0f;
  }
}
