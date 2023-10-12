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

typedef unsigned __int128 uint128_t;

double rlibm_cosf(float x) {
  float_x fX = {.f=x};
  uint32_t b = fX.x<<1;
  if (b < 0xff000000) {
    int k;
    long long a;
    int s = ((fX.x>>23)&0xff) - 150;
    uint64_t m = (fX.x&0x7FFFFF)|1<<23;
    double z, z2;
    if (b < 0x9d000000) {
      if (__builtin_expect(b == 0x9ba56c4e, 0)) return 0x1.9b04988p-7;
      if (__builtin_expect(b == 0x9c27a94a, 0)) return -0x1.8982a08p-26;
      if (b < 0x78921fb6) {
	if (__builtin_expect(b == 0, 0)) return 1.0;
	if (b < 0x73000002) return 0x1.ffffff8p-1;
	z = x;
	z2 = z*z;
	double temp = fma(z2, 0x1.5558e8686159fp-5, -0x1.00000000a7c61p-1);
	return fma(temp, z2, 0x1.0000000000002p+0);
      }
      uint64_t p0 = m*0x441529fc27;
      uint64_t p1 = m*0xa2f9836e4e; p1+=(p0>>40);
      k = (p1>>(33-s));
      a = p1<<(31+s)|((p0<<24)>>(33-s));
    } else {
      if (__builtin_expect(b == 0xaa64a032, 0)) return 0x1.9d4ba48p-1;
      if (__builtin_expect(b == 0xb2887814, 0)) return 0x1.84bec48p-1;
      if (__builtin_expect(b == 0xc22b9622, 0)) return 0x1.f0285d8p-1;
      if (__builtin_expect(b == 0xc2e072ec, 0)) return 0x1.b598ab8p-2;
      if (__builtin_expect(b == 0xe9f4a7f4, 0)) return -0x1.b503d98p-1;
      if (__builtin_expect(b == 0xf2119ae6, 0)) return 0x1.f3176a8p-1;
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
    z = a*0x1p-64;
    z2 = z*z;
    double sinpiK = sinpiMBy256TwoPi[k&511];
    double cospiK = sinpiMBy256TwoPi[(k+128)&511];
    //y=0x1.fffffffffffffp-1 x^(0) + -0x1.3bd3cc9bca02dp-14 x^(2) + 0x1.03c18c953eebep-30 x^(4)
    double cospiZ = fma(z2, 0x1.03c18c953eebep-30, -0x1.3bd3cc9bca02dp-14);
    cospiZ = fma(cospiZ, z2, 0x1.fffffffffffffp-1);
    //y=0x1.921fb54442d39p-7 x^(1) + -0x1.4abbce784a717p-22 x^(3) + 0x1.472122b29fab6p-39 x^(5) 
    double z3 = z2*z;
    z *= 0x1.921fb54442d39p-7;
    double sinpiZ = fma(z2, 0x1.472122b29fab6p-39, -0x1.4abbce784a717p-22);
    sinpiZ = fma(sinpiZ, z3, z);
    return fma(cospiK, cospiZ, -sinpiK*sinpiZ);
  } else { 
    return 0.0f/0.0f;
  }
}
