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
#include <math.h>


#define PI_OVER_512 0x1.921fb54442d18p-8

typedef unsigned __int128 uint128_t;


double rlibm_tanf(float x) { 
  double_x dY;
  float_x fX;
  fX.f = x;
  fX.x &= 0x7FFFFFFF;
  uint64_t sgn = (x < 0) ? 0x8000000000000000 : 0;
  double R, R2;
  if (fX.x < 0x39B89BA3) {
    dY.d = fX.f*0x1.0000008000007p+0;
    dY.x |= sgn;
    return dY.d;
  } else if (fX.x < 0x7F800000) {   
    int N;
    if (fX.x < 0x49800000) {
      if (fX.x < 0x3C490FDB) {
	R = fX.f;
	double R2 = R*R;
	double temp = fma(R2, 0x1.1117b1df3ecb4p-3, 0x1.55555549f55f8p-2);
	dY.d = R*fma(temp, R2, 0x1.0000000000003p+0);
	dY.x |= sgn;
	return dY.d;
      }
      if (fX.x == 0x3F8A1F62) {
	dY.d = 0x1.ddf9f58p+0;
	dY.x ^= sgn;
	return dY.d;
      }
      R = 0x1.45f306cp+7 * fX.f;
      N = R;
      R -= N;
      R = fma(0x1.c9c882a53f84fp-21, fX.f, R);
      if (R > 1.0f) {
	N += 1;
	R -= 1.0f;
      }
    } else {
      if (fX.x == 0x4D56D355) {
	dY.d = 0x1.e803048p-3;
	dY.x ^= sgn;
	return dY.d;
      }
      if (fX.x == 0x5FFD33A4) {
        dY.d = 0x1.a0d9178p+0;
	dY.x ^= sgn;
	return dY.d;
      }
      if (fX.x == 0x6AD36709) {
	dY.d = -0x1.c5612e8p-1;
	dY.x ^= sgn;
	return dY.d;
      }  
      int e = (fX.x >> 23) - 127;
      int s = e - 23;
      uint64_t m = (fX.x&0x7FFFFF)|1<<23;
      static const uint64_t ipi[] = {0xfe5163abdebbc562, 0xdb6295993c439041, 0xfc2757d1f534ddc0, 0xa2f9836e4e441529};
      uint128_t p0 = (uint128_t)m*ipi[0];
      uint128_t p1 = (uint128_t)m*ipi[1]; p1 += p0>>64;
      uint128_t p2 = (uint128_t)m*ipi[2]; p2 += p1>>64;
      uint128_t p3 = (uint128_t)m*ipi[3]; p3 += p2>>64;
      uint64_t p3h = p3>>64, p3l = p3, p2l = p2, p1l = p1;
      uint64_t a;
      if (s < 56) {
	N = (p3h<<(8+s))|(p3l>>(56-s));
	a = (p3l<<(8+s))|(p2l>>(56-s));
      } else if (s == 56) {
	N = p3l;
	a = p2l;
      } else {
	N = (p3l<<(s-56))|(p2l>>(120-s));
	a = (p2l<<(s-56))|(p1l>>(120-s));
      }
      R = a*0x1p-64;
    }
    R *= PI_OVER_512;

    double sinpiM, cospiM, sinpiR, cospiR, temp1, temp2;
    unsigned N2, I;

    double sinR = R;

    N2 = N & 0xFF;
    I = N >> 8;
  	  
    if (I & 0x1) {
      N2 = 255 - N2;
      sinR = PI_OVER_512 - sinR;
    } 
   
    if (I & 0x2) sgn ^= 0x8000000000000000;
   
    double sinR2 = sinR * sinR;

    sinpiM = sinpiMBy512[N2];
    cospiM = cospiMBy512[N2];
    sinpiR, cospiR;

   
    temp1 = fma(sinR2, 0x1.110dd89739b5cp-7, -0x1.55555554ce812p-3);
    sinpiR = sinR*fma(temp1, sinR2, 0x1.ffffffffffffep-1);

    temp2 = fma(sinR2, 0x1.5553dd3610b9p-5, -0x1.ffffffff9c28fp-2);
    cospiR = fma(temp2, sinR2, 0x1.fffffffffffffp-1);


    double sin_x = fma(sinpiM, cospiR, cospiM * sinpiR);

    double cosR = R;

    N2 = N & 0xFF;
    I += 1;
    if (I & 0x2) sgn ^= 0x8000000000000000;

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
 
    temp1 = fma(cosR2, 0x1.110dd89739b5cp-7, -0x1.55555554ce812p-3);
    sinpiR = cosR*fma(cosR2, temp1, 0x1.ffffffffffffep-1);
     
    temp2 = fma(cosR2, 0x1.5553dd3610b9p-5, -0x1.ffffffff9c28fp-2);
    cospiR = fma(temp2, cosR2, 0x1.fffffffffffffp-1);

    double cos_x = fma(sinpiM, sinpiR, cospiM*cospiR);

    dY.d = sin_x/cos_x;
    dY.x ^= sgn;
    return dY.d;
  } else {
    return 0.0f/0.0f;
  }
}
