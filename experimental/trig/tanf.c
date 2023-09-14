#include <stdint.h>
#include <math.h>

#include "pi.h"
#include "rlibm.h"

/* rlibm_tan_pi_256_full_fma_rno */

float tanf(float x) {
  float_x fX = {.f=x};
  uint32_t b = fX.x<<1;
  if (b < 0xff000000) {
    int k;
    double z, z2;
    double sinsinpiK, sincospiK, cossinpiK, coscospiK;
    if (b < 0x96000000) {
      if (b < 0x83000000) {
	if (__builtin_expect(b == 0x7f143ec4, 0)) return 0x1.ddf9f58p+0*__builtin_copysign(1.0f, x);
	if (__builtin_expect(b == 0x8091e348, 0)) return -0x1.e36a458p-10*__builtin_copysign(1.0f, x);
	if (__builtin_expect(b == 0x812d97ac, 0)) return 0x1.2518508p+17*__builtin_copysign(1.0f, x);
	if (__builtin_expect(b == 0x81f6a7a2, 0)) return 0x1.b6e0bf8p+22*__builtin_copysign(1.0f, x);
	if (__builtin_expect(b == 0x8191fef0, 0)) return -0x1.062a398p-9*__builtin_copysign(1.0f, x);
	if (__builtin_expect(b == 0x822d97c8, 0)) return 0x1.99bc5b8p-26*__builtin_copysign(1.0f, x);
	if (__builtin_expect(b == 0x822d98e4, 0)) return 0x1.1c0cce8p-13*__builtin_copysign(1.0f, x);
	if (__builtin_expect(b == 0x822db1be, 0)) return 0x1.9f61278p-9*__builtin_copysign(1.0f, x);
	if (__builtin_expect(b == 0x825fd976, 0)) return 0x1.c068dc8p+11*__builtin_copysign(1.0f, x);
	if (__builtin_expect(b == 0x825fdbd0, 0)) return -0x1.d707988p+16*__builtin_copysign(1.0f, x);
	if (__builtin_expect(b == 0x82c463ac, 0)) return -0x1.aa86798p+24*__builtin_copysign(1.0f, x);
	if (__builtin_expect(b == 0x82d6e936, 0)) return -0x1.87c3508p+0*__builtin_copysign(1.0f, x);
	if (b < 0x78921fb6) {
	  if (b < 0x73713746) return x*0x1.0000008000007p+0;
	  z = x;
	  z2 = z*z;
	  double temp = fma(z2, 0x1.1117b1df3ecb4p-3, 0x1.55555549f55f8p-2);
	  return z*fma(temp, z2, 0x1.0000000000003p+0);
	}
	double idh = 0x1.45f306dc9c883p+6*x, id = __builtin_roundeven(idh);
	double_x dk = {.d = 0x1.8p52 + id};
	k = dk.x;
	z = idh - id;
      } else {
	double dx = x;
	double idl = 0x1.9391054a7f09dp-23*dx, idh = 0x1.45f306dp+6*dx;
	double id = idh + idl;
	double id_int = __builtin_roundeven(id);
	double_x dk = {.d = 0x1.8p52 + id_int};
	k = dk.x;
	z = (idh - id_int) + idl;
      }  
      sinsinpiK = sinpiMBy256TwoPi[k&511];
      sincospiK = sinpiMBy256TwoPi[(k+128)&511];
      cossinpiK = sinsinpiK;
      coscospiK = sincospiK;
    } else {
      if (__builtin_expect(b == 0x9aada6aa, 0)) return 0x1.e803048p-3*__builtin_copysign(1.0f, x);
      if (__builtin_expect(b == 0xafaf61da, 0)) return 0x1.60d1c78p-2*__builtin_copysign(1.0f, x);
      if (__builtin_expect(b == 0xbffa6748, 0)) return 0x1.a0d9178p+0*__builtin_copysign(1.0f, x);
      if (__builtin_expect(b == 0xc7f90dfc, 0)) return 0x1.597f9c8p-1*__builtin_copysign(1.0f, x);
      if (__builtin_expect(b == 0xd5a6ce12, 0)) return -0x1.c5612e8p-1*__builtin_copysign(1.0f, x);
      long long a;
      int s = ((fX.x>>23)&0xff) - 150;
      uint64_t m = (fX.x&0x7FFFFF)|1<<23;
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
      z = a*0x1p-64;
      int sgn = (fX.x>>31)<<8;
      sinsinpiK = sinpiMBy256TwoPi[(k+sgn)&511];
      sincospiK = sinpiMBy256TwoPi[(k+sgn+128)&511];
      cossinpiK = sinpiMBy256TwoPi[k&511];
      coscospiK = sinpiMBy256TwoPi[(k+128)&511];
    }
    z2 = z*z;
    //y=0x1p+0 x^(0) + -0x1.3bd3cc9bf1f2p-14 x^(2) + 0x1.03c20e92da994p-30 x^(4)
    double cospiZ = fma(z2, 0x1.03c20e92da994p-30, -0x1.3bd3cc9bf1f2p-14);
    cospiZ = fma(cospiZ, z2, 0x1p+0);
    //y=0x1.921fb54442d1cp-7 x^(1) + -0x1.4abbce65b76d9p-22 x^(3) + 0x1.46911d1f4347cp-39 x^(5)
    double z3 = z2*z;
    z *= 0x1.921fb54442d1cp-7;
    double sinpiZ = fma(z2, 0x1.46911d1f4347cp-39, -0x1.4abbce65b76d9p-22);
    sinpiZ = fma(sinpiZ, z3, z);
    return fma(sinsinpiK, cospiZ, sincospiK*sinpiZ)/fma(coscospiK, cospiZ, -cossinpiK*sinpiZ);
  } else { 
    return 0.0f/0.0f;
  }
}
