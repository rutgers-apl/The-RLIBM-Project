/*
MIT License

Copyright (c) 2022 Sehyeok Park, Mridul Aanjaneya, and Santosh Nagarakatte, Rutgers
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



#include "rlibm.h"
#include "sinpicospi.h"
#include "pi.h"
#include "math.h"

#define PI_OVER_256 0x1.921fb54442d18p-7 

typedef unsigned __int128 uint128_t;

double rlibm_cosf(float x) {
  double_x dY;
  float_x fX;
  fX.f = x;
  fX.x &= 0x7FFFFFFF;
  uint64_t sgn = 0;
  double R, R2;
  if (fX.x == 0) {
    return 1.0;
  } 
  if (fX.x < 0x39800001) {
    return 0x1.fffffffffffffp-1;
  } else if (fX.x < 0x7F800000) {
    int N;
    if (fX.x < 0x4A000000) {
      if (fX.x < 0x3C490FDB) {
	R = fX.f;
	R2 = R*R;
	double temp = fma(R2, 0x1.5558e8686159fp-5, -0x1.00000000a7c61p-1);
	return fma(temp, R2, 0x1.0000000000002p+0);
      }
      R = 0x1.45f306cp+6 * fX.f;
      N = R;
      R -= N;
      R = fma(0x1.c9c882a53f84fp-22, fX.f, R);
      if (R > 1.0f) {
	N += 1;
	R -= 1.0f;
      }
    } else {
      if (fX.x == 0x5922AA80) {
	return 0x1.115d7d8p-1;
      }
      if (fX.x == 0x6A127977) {
	return 0x1.af5c698p-2;
      }
      if (fX.x == 0x7908CD73) {
	return 0x1.f3176a8p-1;
      }
      if (fX.x == 0x7A38AB34) {
	return 0x1.f663298p-1;
      }
      int e = (fX.x >> 23) - 127;
      int s = e - 23;
      uint64_t m = (fX.x&0x7FFFFF)|1<<23;

      /* 256 bits of 1/pi in 4 64-bit integers */
      static const uint64_t ipi[] = {0xfe5163abdebbc562, 0xdb6295993c439041, 0xfc2757d1f534ddc0, 0xa2f9836e4e441529};
      uint128_t p0 = (uint128_t)m*ipi[0];
      uint128_t p1 = (uint128_t)m*ipi[1]; p1 += p0>>64;
      uint128_t p2 = (uint128_t)m*ipi[2]; p2 += p1>>64;
      uint128_t p3 = (uint128_t)m*ipi[3]; p3 += p2>>64;
      uint64_t p3h = p3>>64, p3l = p3, p2l = p2, p1l = p1;
      uint64_t a;
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
      R = a*0x1p-64;
    }
    R *= PI_OVER_256;

    double sinpiM;
    double cospiM;
    
    unsigned N2 = N & 0x7F;
    unsigned I = (N >> 7) + 1;
    if (I & 0x2) sgn ^= 0x8000000000000000;

    if (I & 1) {
      if (N2 == 0) {
	cospiM = 1.0;
	sinpiM = 0.0;
      } else {
	N2++;
	R = PI_OVER_256 - R;
	cospiM = sinpiMBy256[128 - N2];
	sinpiM = cospiMBy256[128 - N2];
      }
    } else {
      cospiM = sinpiMBy256[N2];
      sinpiM = cospiMBy256[N2];
    }
    R2 = R*R;
    double sinpiR, cospiR, temp1, temp2; 

    temp1 = fma(R2, 0x1.1111de524b6fp-7, -0x1.55555555d376p-3);
    sinpiR = R*fma(R2, temp1, 0x1.0000000000001p+0);
     
    temp2 = fma(R2, 0x1.555488594da9dp-5, -0x1.ffffffff83643p-2);
    cospiR = fma(R2, temp2, 0x1.ffffffffffffcp-1);
    
    dY.d = fma(sinpiM, sinpiR, cospiM*cospiR);
    dY.x ^= sgn;
    return dY.d;
  } else {
    return 0.0f/0.0f;
  }
}
