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
#include "pi.h"


#define PI_OVER_256 0x1.921fb54442d18p-7

typedef unsigned __int128 uint128_t;

double rlibm_sinf(float x) {
  float_x fX = {.f=x};
  uint32_t b = fX.x<<1;
  if (b < 0xff000000) {
    int k;
    long long a;
    int s = ((fX.x>>23)&0xff) - 150;
    uint64_t m = (fX.x&0x7FFFFF)|1<<23;
    double z, z2;
    if (b < 0x9d000000) {
      if (__builtin_expect(b == 0x95fbd9c8, 0)) {
	double y = -0x1.ff6dc18p-1*__builtin_copysign(1.0, x);
	return y;
      }
      if (__builtin_expect(b == 0x995c4a68, 0)) {
	double y = -0x1.181a858p-6*__builtin_copysign(1.0, x);
	return y;
      }
      if (__builtin_expect(b == 0x9cb0e51c, 0)) {
	double y = -0x1.851fd78p-14*__builtin_copysign(1.0, x);
	return y;
      }
      if (b < 0x79921fb6) {
	if (b < 0x73d12ed2) return x*0x1.ffffffffffffep-1;
	z = x;
	z2 = z*z;
	double temp = fma(z2, 0x1.11130658ac4e4p-7, -0x1.55555558c07f5p-3);
	return z*fma(temp, z2, 0x1.0000000000011p+0);
      }
      uint64_t p0 = m*0x441529fc27;
      uint64_t p1 = m*0xa2f9836e4e; p1+=(p0>>40);
      k = (p1>>(33-s));
      a = p1<<(31+s)|((p0<<24)>>(33-s));
    } else {
      if (__builtin_expect(b == 0xa4728fec, 0)) {
	double y = -0x1.24f23b8p-1*__builtin_copysign(1.0, x);
	return y;
      }
      if (__builtin_expect(b == 0xab95f654, 0)) {
	double y = -0x1.fcf42d8p-1*__builtin_copysign(1.0, x);
	return y;
      }
      if (__builtin_expect(b == 0xcf524856, 0)) {
	double y = -0x1.ff57018p-1*__builtin_copysign(1.0, x);
	return y;
      }
      if (__builtin_expect(b == 0xe6487e0c, 0)) {
	double y = 0x1.2875078p-2*__builtin_copysign(1.0, x);
	return y;
      }
      static const uint64_t ipi[] = {0xdb6295993c439041, 0xfc2757d1f534ddc0, 0xa2f9836e4e441529};
      uint128_t p1 = (uint128_t)m*ipi[0];
      uint128_t p2 = (uint128_t)m*ipi[1]; p2 += p1>>64;
      uint128_t p3 = (uint128_t)m*ipi[2]; p3 += p2>>64;
      uint64_t p3h = p3>>64, p3l = p3, p2l = p2, p1l = p1;
      if (s < 57) {
	k = (p3h<<(7+s))|(p3l>>(57-s));
	a = (p3l<<(7+s))|(p2l>>(57-s));
      } else if (s == 57) {
	k = p3l;
	a = p2l;
      } else {
	k = (p3l<<(s-57))|(p2l>>(121-s));
	a = (p2l<<(s-57))|(p1l>>(121-s));
      }
    }
    long sm = a>>63;
    k -= sm;
    k += (fX.x>>31)<<8;
    z = a*0x1p-64;
    z2 = z*z;
    double sinpiK = sinpiMBy256TwoPi[k&511];
    double cospiK = sinpiMBy256TwoPi[(k+128)&511];
    //y=0x1.fffffffffffffp-1 x^(0) + -0x1.3bd3cc9baf37fp-14 x^(2) + 0x1.03c102099750dp-30 x^(4)
    double cospiZ = fma(z2, 0x1.03c102099750dp-30, -0x1.3bd3cc9baf37fp-14);
    cospiZ = fma(cospiZ, z2, 0x1.fffffffffffffp-1);
    double z3 = z2*z;
    //y=0x1.921fb54442d1cp-7 x^(1) + -0x1.4abbce67633e5p-22 x^(3) + 0x1.46a562b83e363p-39 x^(5)
    z *= 0x1.921fb54442d1cp-7;
    double sinpiZ = fma(z2, 0x1.46a562b83e363p-39, -0x1.4abbce67633e5p-22);
    sinpiZ = fma(sinpiZ, z3, z);
    return fma(sinpiK, cospiZ, cospiK*sinpiZ);
  } else { 
    return 0.0f/0.0f;
  }
}



