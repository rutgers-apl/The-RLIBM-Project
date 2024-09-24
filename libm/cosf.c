/*
MIT License
Copyright (c) 2024 Sehyeok Park, Mridul Aanjaneya, and Santosh
Nagarakatte, Rutgers Architecture and Programming Languages (RAPL)
Group

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions: The above copyright notice and this
permission notice shall be included in all copies or substantial
portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/


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
      if (b < 0x78921fb6) {
	if (__builtin_expect(b == 0, 0)) return 1.0;
	if (b < 0x73000002) return 0x1.ffffff8p-1;
	//y=0x1.0000000000002p+0 x^(0) + -0x1.00000000752fdp-1 x^(2) + 0x1.5557b658a7321p-5 x^(4)
	z = x;
	z2 = z*z;
	double temp = __builtin_fma(z2, 0x1.5557b658a7321p-5, -0x1.00000000752fdp-1);
	return __builtin_fma(z2, temp, 0x1.0000000000002p+0);
      }
      uint64_t p0 = m*0x441529fc27;
      uint64_t p1 = m*0xa2f9836e4e; p1+=(p0>>40);
      k = (p1>>(33-s));
      a = p1<<(31+s);
      if (b > 0x8b400000) a |= ((p0<<24)>>(33-s));
      long sm = a>>63;
      k -= sm;
      z = ((a>>10)<<10)*0x1p-64;
      z2 = z*z;
      double sinpiK = sinpiMBy256TwoPi[k&511];
      double cospiK = sinpiMBy256TwoPi[(k+128)&511];
      double cospiZ;
      if (__builtin_expect(z == -0x1.6ad140905564p-6, 0)) cospiZ = 0x1.fffffe2611ef2p-1;
      else if (__builtin_expect(z == 0x1.868f3be09e38p-2, 0)) cospiZ = 0x1.fffe91699f62cp-1;
      else {
	//-0x1.6ad140905564p-6, 0x1.868f3be09e38p-2
	//y=0x1.fffffffffffffp-1 x^(0) + -0x1.3bd3cc9c623a3p-14 x^(2) + 0x1.03c42bb97e0aep-30 x^(4) 
	cospiZ = __builtin_fma(z2, 0x1.03c42bb97e0aep-30, -0x1.3bd3cc9c623a3p-14);
	cospiZ = __builtin_fma(cospiZ, z2, 0x1.fffffffffffffp-1);
      }
      double sinpiZ;
      if (__builtin_expect(z == 0x1.f50866898p-20, 0)) sinpiZ = 0x1.8982a08p-26;
      else if (__builtin_expect(z == 0x1.d99d2c356a4p-5, 0)) sinpiZ = 0x1.73f9baeba657fp-11;
      else {
	//0x1.f50866898p-20, 0x1.d99d2c356a4p-5
	//y=0x1.921fb54442d06p-7 x^(1) + -0x1.4abbce622bdf3p-22 x^(3) + 0x1.469f84cd0a117p-39 x^(5)
	sinpiZ = __builtin_fma(z2, 0x1.469f84cd0a117p-39, -0x1.4abbce622bdf3p-22);
	sinpiZ = z*(__builtin_fma(sinpiZ, z2, 0x1.921fb54442d06p-7));
      }
      return __builtin_fma(cospiK, cospiZ, -sinpiK*sinpiZ);
    } else {
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
      long sm = a>>63;
      k -= sm;
      z = ((a>>10)<<10)*0x1p-64;
      z2 = z*z;
      double sinpiK = sinpiMBy256TwoPi[k&511];
      double cospiK = sinpiMBy256TwoPi[(k+128)&511];
      double cospiZ;
      if (__builtin_expect(z == -0x1.bb0c2d47cfd2cp-2, 0)) cospiZ = 0x1.fffe266eed5ebp-1;
      else if (__builtin_expect(z == -0x1.5d823ac42b395p-2, 0)) cospiZ = 0x1.fffed9124fb54p-1;
      else if (__builtin_expect(z == -0x1.0788e4762434p-2, 0)) cospiZ = 0x1.ffff57fef756ep-1;
      else if (__builtin_expect(z == -0x1.a9e8818a3af1ap-3, 0)) cospiZ = 0x1.ffff91c2bf126p-1;
      else if (__builtin_expect(z == -0x1.aafc551c92bf8p-4, 0)) cospiZ = 0x1.ffffe579a26fp-1;
      else if (__builtin_expect(z == -0x1.7dcb0e32af9b4p-4, 0)) cospiZ = 0x1.ffffe98bcf4c4p-1;
      else if (__builtin_expect(z == -0x1.3d5cf3c8209fp-6, 0)) cospiZ = 0x1.ffffffa2ed8e4p-1;
      else if (__builtin_expect(z == 0x1.8155d60fb924p-3, 0)) cospiZ = 0x1.ffffa60c321cp-1;
      else if (__builtin_expect(z == 0x1.58f335ff0f9d9p-2, 0)) cospiZ = 0x1.fffee0c5cde96p-1;
      else {	  
	//-0x1.bb0c2d47cfd2cp-2, -0x1.5d823ac42b395p-2, -0x1.0788e4762434p-2, -0x1.a9e8818a3af1ap-3, -0x1.aafc551c92bf8p-4, -0x1.7dcb0e32af9b4p-4, -0x1.3d5cf3c8209fp-6, 0x1.8155d60fb924p-3, 0x1.58f335ff0f9d9p-2
	//y=0x1.fffffffffffffp-1 x^(0) + -0x1.3bd3cc9bc5f2ap-14 x^(2) + 0x1.03c24f83d10aap-30 x^(4)
	cospiZ = __builtin_fma(z2, 0x1.03c24f83d10aap-30, -0x1.3bd3cc9bc5f2ap-14);
	cospiZ = __builtin_fma(cospiZ,z2,0x1.fffffffffffffp-1);
      }
      //y=0x1.921fb54442d3ap-7 x^(1) + -0x1.4abbce79b84cp-22 x^(3) + 0x1.472da2d3e540ep-39 x^(5)
      double sinpiZ = __builtin_fma(z2, 0x1.472da2d3e540ep-39, -0x1.4abbce79b84cp-22);
      sinpiZ = z*__builtin_fma(sinpiZ,z2,0x1.921fb54442d3ap-7);
      return __builtin_fma(cospiK,cospiZ,-sinpiK*sinpiZ);
    }
  } else { 
    return 0.0f/0.0f;
  }
}

float cr_cosf(float x){
  return (float) rlibm_cosf(x);
}

float cosf(float x){
  return (float) rlibm_cosf(x);
}
