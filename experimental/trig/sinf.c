#include <stdint.h>
#include <math.h>

#include "pi.h"
#include "rlibm.h"

/* rlibm_sin_pi_256_full_fma_rno */


float sinf(float x) {
  float_x fX = {.f=x};
  uint32_t b = fX.x<<1;
  if (b < 0xff000000) {
    int k;
    double z, z2;
    if (b < 0x96000000) {
      if (b < 0x83000000) {
	if (b < 0x79921fb6) {
	  if (b < 0x73d12ed2) return x*0x1.ffffffffffffep-1;
	  z = x;
	  z2 = z*z;
	  double temp = fma(z2, 0x1.11130658ac4e4p-7, -0x1.55555558c07f5p-3);
	  return z*fma(temp, z2, 0x1.0000000000011p+0);
	}
	if (__builtin_expect(b == 0x7e75b8a2, 0)) return 0x1.5568898p-1*__builtin_copysign(1.0f, x);
	if (__builtin_expect(b == 0x7f4f0654, 0)) return 0x1.ee836b8p-1*__builtin_copysign(1.0f, x);
	if (__builtin_expect(b == 0x812d97c8, 0)) {
	  double y = -0x1.ffffff8p-1*__builtin_copysign(1.0f, x);
	  return y;
	}
	if (__builtin_expect(b == 0x822d97c8, 0)) return -0x1.99bc5b8p-26*__builtin_copysign(1.0f, x);
	if (__builtin_expect(b == 0x82a41896, 0)) return 0x1.10aca08p-1*__builtin_copysign(1.0f, x);
	if (__builtin_expect(b == 0x82f66f3c, 0)) return 0x1.c333c18p-8*__builtin_copysign(1.0f, x);
	double idh = 0x1.45f306dc9c883p+6*x, id = __builtin_roundeven(idh);
	double_x dk = {.d = 0x1.8p52 + id};
	k = dk.x;
	z = idh - id;
      } else {
        if (__builtin_expect(b == 0x86f9cbe2, 0)) {
	  double y = 0x1.ffffff8p-1*__builtin_copysign(1.0f, x);
	  return y;
	}
	if (__builtin_expect(b == 0x95fbd9c8, 0)) return -0x1.ff6dc18p-1*__builtin_copysign(1.0f, x);
	double dx = x;
	double idl = 0x1.9391054a7f09dp-23*dx, idh = 0x1.45f306dp+6*dx;
	double id = idh + idl;
	double id_int = __builtin_roundeven(id);
	double_x dk = {.d = 0x1.8p52 + id_int};
	k = dk.x;
	z = (idh - id_int) + idl;
      }
    } else {
      if (__builtin_expect(b == 0x984665d2, 0)) {
        double y = -0x1.ffffff8p-1*__builtin_copysign(1.0f, x);
	return y;
      }
      if (__builtin_expect(b == 0xa147d0fe, 0)) {
	double y = 0x1.ffffff8p-1*__builtin_copysign(1.0f, x);
	return y;
      }
      if (__builtin_expect(b == 0xa4728fec, 0)) return -0x1.24f23b8p-1*__builtin_copysign(1.0f, x);
      if (__builtin_expect(b == 0xa7628d4c, 0)) {
	double y = -0x1.ffffff8p-1*__builtin_copysign(1.0f, x);
	return y;
      }
      if (__builtin_expect(b == 0xb80f79a0, 0)) {
	double y = -0x1.ffffff8p-1*__builtin_copysign(1.0f, x);
	return y;
      }
      if (__builtin_expect(b == 0xcb130930, 0)) {
	double y = 0x1.ffffff8p-1*__builtin_copysign(1.0f, x);
	return y;
      }
      if (__builtin_expect(b == 0xcf524856, 0)) return -0x1.ff57018p-1*__builtin_copysign(1.0f, x);
      if (__builtin_expect(b == 0xd432ede2, 0)) {
	double y = 0x1.ffffff8p-1*__builtin_copysign(1.0f, x);
	return y;
      }
      if (__builtin_expect(b == 0xd8abb4b0, 0)) {
	double y = 0x1.ffffff8p-1*__builtin_copysign(1.0f, x);
	return y;
      }
      if (__builtin_expect(b == 0xdef37c8a, 0)) {
	double y = 0x1.ffffff8p-1*__builtin_copysign(1.0f, x);
	return y;
      }
      if (__builtin_expect(b == 0xeeb08c4a, 0)) {
	double y = -0x1.ffffff8p-1*__builtin_copysign(1.0f, x);
	return y;
      }
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
      k += (fX.x>>31)<<8;
      z = a*0x1p-64;
    }
    z2 = z*z;
    double sinpiK = sinpiMBy256TwoPi[k&511];
    double cospiK = sinpiMBy256TwoPi[(k+128)&511];
    //y=0x1p+0 x^(0) + -0x1.3bd3cc9be1f91p-14 x^(2) + 0x1.03c1e4ff8a2d1p-30 x^(4)
    double cospiZ = fma(z2, 0x1.03c1e4ff8a2d1p-30, -0x1.3bd3cc9be1f91p-14);
    cospiZ = fma(cospiZ, z2, 0x1p+0);
    double z3 = z2*z;
    //y=0x1.921fb54442d1cp-7 x^(1) + -0x1.4abbce66d3629p-22 x^(3) + 0x1.469fe73456b22p-39 x^(5)
    z *= 0x1.921fb54442d1cp-7;
    double sinpiZ = fma(z2, 0x1.469fe73456b22p-39, -0x1.4abbce66d3629p-22);
    sinpiZ = fma(sinpiZ, z3, z);
    return fma(sinpiK, cospiZ, cospiK*sinpiZ);
  } else { 
    return 0.0f/0.0f;
  }
}
