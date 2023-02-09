/*
MIT License
Copyright (c) 2023 Sehyeok Park, Mridul Aanjaneya, and Santosh
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


double rlibm_cosf(float x) {  
  double y;
  float_x fX;
  fX.f = x;
  fX.x &= 0x7FFFFFFF;
  if (fX.x == 0) {
    return 1.0;
  }
  if (fX.x < 0x39800001) {
    y = 0x1.fffffffffffffp-1;
    return y;
  } else if (fX.x < 0x7F800000) {
    double R, R2;
    if (fX.x < 0x3C490FDB) {
      R = fX.f;
      R2 = R*R;
      double temp = fma(R2, 0x1.5558e8686159fp-5, -0x1.00000000a7c61p-1);
      y = fma(temp, R2, 0x1.0000000000002p+0);
      return y;
    }
    int s = (fX.x >> 23) - 150;
    uint64_t m = (fX.x&0x7FFFFF)|1<<23;
    int N;
    uint64_t a;    
    if (fX.x < 0x4A000000) {
      if (fX.x == 0x424790ce) {
	y = 0x1.dc98028p-1;
	return y;
      }
      uint64_t p0 = m*0x441529fc27;
      uint64_t p1 = m*0xa2f9836e4e; p1+=(p0>>40);
      N = (p1>>(33-s));
      a = p1<<(31+s)|((p0<<24)>>(33-s));
    } else {
      if (fX.x == 0x6A127977) {
	y = 0x1.af5c698p-2;
	return y;
      }
      if (fX.x == 0x7908CD73) {
	y = 0x1.f3176a8p-1;
	return y;
      }
      if (fX.x == 0x7A38AB34) {
	y = 0x1.f663298p-1;
	return y;
      }
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
    double sinpiM, cospiM;
    long N2 = N & 0x7F;
    unsigned I = (N >> 7) + 1;
    if ((I & 1) & (N2 == 0)) {
      cospiM = 1.0;
      sinpiM = 0.0;
    } else {
      a ^= ~(0xFFFFFFFFFFFFFFFF+(I & 0x1));
      N2 ^= ~(0xFFFFFFFF+(I & 0x1));
      N2 += (128&((I & 0x1)<<7));
      cospiM = sinpiMBy256[(unsigned)N2];
      sinpiM = cospiMBy256[(unsigned)N2];
    }
    double temp1, temp2;
    R = a*0x1p-64*PI_OVER_256;
    R2 = R*R;

    temp1 = fma(R2, 0x1.1111de524b6fp-7, -0x1.55555555d376p-3);
    y = sinpiM*R*fma(R2, temp1, 0x1.0000000000001p+0); 
    
    temp2 = fma(R2, 0x1.55546cdabdep-5, -0x1.ffffffff5c8f6p-2);
    y += cospiM*fma(R2, temp2, 0x1.ffffffffffffbp-1);

    double sgns[2] = {1.0, -1.0};
    return y*sgns[(I & 0x2)>>1];
  } else {
    return 0.0f/0.0f;
  }
}
