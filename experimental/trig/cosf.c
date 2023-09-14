#include <stdint.h>
#include <math.h>

#include "pi.h"
#include "rlibm.h"

/* rlibm_cos_pi_256_full_fma_rno */

float cosf(float x) {
  float_x fX = {.f=x};
  uint32_t b = fX.x<<1;
  if (b < 0xff000000) {
    int k;
    double z, z2;
    if (b < 0x96000000) {
      if (b < 0x83000000) {
	if (b < 0x78921fb6) {
          if (__builtin_expect(b == 0, 0)) return 1.0f;
	  if (b < 0x73000002) return 1.0f - 0x1p-25f;
	  z = x;
	  z2 = z*z;
	  double temp = fma(z2, 0x1.5558e8686159fp-5, -0x1.00000000a7c61p-1);
	  return fma(temp, z2, 0x1.0000000000002p+0);
	}
	if (__builtin_expect(b == 0x812d97c8, 0)) return 0x1.99bc5cp-27f - 0x1p-52f;
	if (__builtin_expect(b == 0x82c463ac, 0)) return -0x1.334d44p-25f - 0x1p-50f; 
	double idh = 0x1.45f306dc9c883p+6*x, id = __builtin_roundeven(idh);
	double_x dk = {.d = 0x1.8p52 + id};
	k = dk.x;
	z = idh - id;
      } else {
	if (__builtin_expect(b == 0x87f8cbe2, 0)) return -0x1.14a282p-1f + 0x1p-26f;
	if (__builtin_expect(b == 0x87f9cbe2, 0)) return -1.0f + 0x1p-25f;
	double dx = x;
	double idl = 0x1.9391054a7f09dp-23*dx, idh = 0x1.45f306dp+6*dx;
	double id = idh + idl;
	double id_int = __builtin_roundeven(id);
	double_x dk = {.d = 0x1.8p52 + id_int};
	k = dk.x;
	z = (idh - id_int) + idl;
      }
    } else {
      if (__builtin_expect(b == 0xa247d0fe, 0)) return -1.0f + 0x1p-25f;
      if (__builtin_expect(b == 0xa347d0fe, 0)) return 1.0f - 0x1p-25f;
      if (__builtin_expect(b == 0xa8628d4c, 0)) return -1.0f + 0x1p-25f; 
      if (__builtin_expect(b == 0xaa64a032, 0)) return 0x1.9d4ba4p-1f + 0x1p-26f;
      if (__builtin_expect(b == 0xb2887814, 0)) return 0x1.84bec4p-1f + 0x1p-26f;
      if (__builtin_expect(b == 0xb548a858, 0)) return 0x1.f48148p-2f + 0x1p-27f;
      if (__builtin_expect(b == 0xc22b9622, 0)) return 0x1.f0285ep-1f - 0x1p-26f;
      if (__builtin_expect(b == 0xcc130930, 0)) return -1.0f + 0x1p-25f;
      if (__builtin_expect(b == 0xd532ede2, 0)) return -1.0f + 0x1p-25f;
      if (__builtin_expect(b == 0xdff37c8a, 0)) return -1.0f + 0x1p-25f;
      if (__builtin_expect(b == 0xe0f37c8a, 0)) return 1.0f - 0x1p-25f;
      if (__builtin_expect(b == 0xe9f4a7f4, 0)) return -0x1.b503dap-1f + 0x1p-26f;
      if (__builtin_expect(b == 0xf2119ae6, 0)) return 0x1.f3176ap-1f + 0x1p-26f;
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
    }
    z2 = z*z;
    double sinpiK = sinpiMBy256TwoPi[k&511];
    double cospiK = sinpiMBy256TwoPi[(k+128)&511];
    //y=0x1p+0 x^(0) + -0x1.3bd3cc9be344p-14 x^(2) + 0x1.03c1cf69ce864p-30 x^(4)    
    double cospiZ = fma(z2, 0x1.03c1cf69ce864p-30, -0x1.3bd3cc9be344p-14);
    cospiZ = fma(cospiZ, z2, 0x1p+0);
    //y=0x1.921fb54442d3cp-7 x^(1) + -0x1.4abbce7c857bfp-22 x^(3) + 0x1.474adac8433a9p-39 x^(5) 
    double z3 = z2*z;
    z *= 0x1.921fb54442d3cp-7;
    double sinpiZ = fma(z2, 0x1.474adac8433a9p-39, -0x1.4abbce7c857bfp-22);
    sinpiZ = fma(sinpiZ, z3, z);
    return fma(cospiK, cospiZ, -sinpiK*sinpiZ);
  } else { 
    return 0.0f/0.0f;
  }
}
