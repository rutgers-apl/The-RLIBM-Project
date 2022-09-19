#include <math.h>
#include "rlibm.h"
#include "exp2.h"

#define MAXVAL 3.40282361850336062550457001444955389952e+38
#define MAXm1VAL 3.40282356779733661637539395458142568448e+38

double rlibm_exp10(float x) {
  float_x fx;
  fx.f = x;
  
  // Take care of special cases
  if ((fx.x & 0x7FFFFFFF) == 0) return 1.0;

  if (fx.x <= 861821911) {
    if (fx.x <= 853433304) return 1.0000000298023223876953125;
    return 1.0000000894069671630859375;
  }

  if (1109008539 <= fx.x && fx.x <= 3000916953) {
    if (fx.x < 0x80000000) {
      if (fx.x < 0x7F800000) return 3.40282361850336062550457001444955389952e+38;
      if (fx.x == 0x7F800000) return 1.0 / 0.0;
      return 0.0/0.0;
    }

    // negative counterpart
    if (fx.x <= 2992528344) return 0.99999998509883880615234375;

    return 0.99999995529651641845703125;
  }

  if (fx.x >= 3258228278) {
    if (fx.x == 0xFF800000) return 0.0;
    if (fx.x < 0xFF800000) return ldexp(1.0, -151);
    return 0.0/0.0;
  }

  // If x == 0.0, 1.0, 2.0, ..., 10.0, then it's also special case
  switch(fx.x) {
  case 0x00000000:
  case 0x80000000: return 1.0;
  case 0x3f800000: return 10.0;
  case 0x40000000: return 100.0;
  case 0x40400000: return 1000.0;
  case 0x40800000: return 10000.0;
  case 0x40a00000: return 100000.0;
  case 0x40c00000: return 1000000.0;
  case 0x40e00000: return 10000000.0;
  case 0x41000000: return 100000000.0;
  case 0x41100000: return 1000000000.0;
  case 0x41200000: return 10000000000.0;
  }
  
  // Perform range reduction
  double xp = x * 2.12603398072791179629348334856331348419189453125e+02;
  int N = (int)xp;
  int N2 = N % 64;
  if (N2 < 0) N2 += 64;
  int N1 = N - N2;
  
  int M = N1 / 64;
  int J = N2;
  double R = x - N *
  4.703593682249706219022922226713490090332925319671630859375e-03;
  
  double_x dx;
  dx.d = R;
  // Find the polynomial coefficients to use.
  double y;

  if(R == -0x1.0c34ef4921p-9){
    y = 0x1.fd97e23938223p-1;
  }else if (R == -0x1.82c1c01d78p-10){
    y = 0x1.fe437ac045045p-1;
  }else if(R == 0x1.fb2943886481p-10){
    y = 0x1.012498c03e2ep+0;
  }else if(R == 0x1.d924c938a8p-9){
    y = 0x1.0222fe9de751bp+0;
  }else{

    double coeffs[]={
		     0x1p+0,
		     0x1.26bb1bbb5575p+1,
		     0x1.53524c7378f8cp+1,
		     0x1.04705809a324ap+1,
		     0x1.2bd81086fea0fp+0,
		     0x1.182ea56fde13cp-1
    };

    double xsquare = R * R;
    double temp1 = fma(R, coeffs[1], coeffs[0]);
    double temp2 = fma(R, coeffs[5], coeffs[4]);
    double temp3 = fma(R, coeffs[3], coeffs[2]);

    double temp4 = fma(xsquare, temp2, temp3);
    y= fma(xsquare, temp4, temp1);

  }
  
  // Perform output compensation
  y *= ldexp(exp2JBy64[J], M);
  return y;
}
