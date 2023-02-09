#include <math.h>
#include "rlibm.h"
#include "sinpicospi.h"
#include "pi.h"

typedef unsigned __int128 uint128_t;


double rlibm_tanf(float x) {  
  float_x fX;
  fX.f = x;
  int64_t sb = fX.x>>31;
  double sgns[2] = {1.0, -1.0};
  double sgn = sgns[sb];
  fX.x &= 0x7FFFFFFF;
  double y;
  if (fX.x < 0x39B89BA3) {
    y = fX.f*0x1.0000008000007p+0;
    return y*sgn;
  } else if (fX.x < 0x7F800000) {
    double R, R2;
    if (fX.x < 0x3C490FDB) {
      R = fX.f;
      R2 = R*R;
      double temp = fma(R2, 0x1.1117b1df3ecb4p-3, 0x1.55555549f55f8p-2);
      y = R*fma(temp, R2, 0x1.0000000000003p+0);
      return y*sgn;
    }
    int e = (fX.x >> 23) - 127;
    int s = e - 23;
    uint64_t m = (fX.x&0x7FFFFF)|1<<23;
    uint64_t N;
    int64_t a;
    if (fX.x < 0x4D800000) {
      if (fX.x == 0x3f2aed4f) {
	return 0x1.93b50c8p-1*sgn;
      }
      if (fX.x == 0x408174dd) {
	return 0x1.4536628p+0*sgn;
      }
      if (fX.x == 0x40e67f59) {
	return 0x1.5019108p+0*sgn;
      }
      if (fX.x == 0x4d16f93a) {
	return -0x1.dd4a7d8p-16*sgn;
      }
      if (fX.x == 0x4d56d355) {
	return 0x1.e803048p-3*sgn;
      }
      uint64_t p0 = m*0x441529fc27;
      uint64_t p1 = m*0xa2f9836e4e; p1+=(p0>>40);
      uint64_t g = (0xffffffffffffffff + (s>=-23));
      a = (g & (p1>>(-24-s))) | (~g & (((p0<<24)>>(40-s))|(p1<<(24+s))));
      N = (~g & (p1>>(40-s)));
    } else {
      if (fX.x == 0x5980445e) {
        return 0x1.ca1ed08p+0*sgn;
      }
      if (fX.x == 0x5ffd33a4) {
	return 0x1.a0d9178p+0*sgn;
      }
      if (fX.x == 0x63fc86fe) {
	return 0x1.597f9c8p-1*sgn;
      }
      if (fX.x == 0x72b505bb) {
	return -0x1.e42a1e8p+0*sgn;
      }
      if (fX.x == 0x77cda26b) {
	return 0x1.1056678p+0*sgn;
      }
      if (fX.x == 0x79c42c65) {
	return 0x1.45c66c8p+0*sgn;
      }
      if (fX.x == 0x7bedc1ea) {
	return 0x1.403ed58p+0*sgn;
      }
      if (fX.x == 0x7c8d49d9) {
	return 0x1.06d84e8p+0*sgn;
      }
      static const uint64_t ipi[] = {0xfe5163abdebbc562, 0xdb6295993c439041, 0xfc2757d1f534ddc0, 0xa2f9836e4e441529};
      uint128_t p0 = (uint128_t)m*ipi[0];
      uint128_t p1 = (uint128_t)m*ipi[1]; p1 += p0>>64;
      uint128_t p2 = (uint128_t)m*ipi[2]; p2 += p1>>64;
      uint128_t p3 = (uint128_t)m*ipi[3]; p3 += p2>>64;
      uint64_t p3h = p3>>64, p3l = p3, p2l = p2, p1l = p1;
      if(s<64) {
	N = (p3h<<s)|(p3l>>(64-s));
	a = p3l<<s|p2l>>(64-s);
      } else if(s==64) {
	N = p3l;
	a = p2l;
      } else {
	N = p3l<<s|p2l>>(128-s);
	a = p2l<<s|p1l>>(128-s);
      }
    }
    N += ((uint64_t)a>>63);
    sb = (sb<<63)>>63;
    N = (N^sb) - sb;
    R = (a^sb)*0x1p-64;
    R2 = R*R;
    double R4 = R2*R2, R8 = R4*R4, R9 = R8*R;

    double temp1 = fma(-0x1.4abbce625bb1dp-1, R2, 0x1.921fb54442d16p+0);
    double temp2 = fma(-0x1.32d2cc8ac7a3cp-8, R2, 0x1.466bc6770bfebp-4);  
    double temp3 = fma(-0x1.e2aba788c8129p-19, R2, 0x1.5077d5149c381p-13);  
    temp1 = fma(temp2, R4, temp1);
    temp3 = fma(0x1.c6a694fc3b703p-25, R4, temp3); 
    double sinR = fma(temp3, R9, temp1*R);

    double temp4 = fma(-0x1.3bd3cc9be45d1p+0, R2, 0x1.fffffffffffffp-1);
    double temp5 = fma(-0x1.55d3c7e41fb1bp-6, R2, 0x1.03c1f081b6166p-2); 
    double temp6 = fma(-0x1.a6ceddbee43ebp-16, R2, 0x1.e1f5051c49684p-11); 
    temp4 = fma(temp5, R4, temp4);
    temp6 = fma(0x1.f5cdb3ae720ap-22, R4, temp6); 
    double cosR = fma(temp6, R8, temp4);
    
    if (N%2) {
      return -cosR/sinR;
    } else {
      return sinR/cosR;
    }
  } else {
    return 0.0f/0.0f;
  }
}
