#include <math.h>
#include "float-rno/float_rno_lib.h"
#include "exp2.h"

#define MAXVAL 3.40282361850336062550457001444955389952e+38
#define MAXm1VAL 3.40282356779733661637539395458142568448e+38

double rlibm_rno_exp10(float x) {
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
  
  switch (dx.x) {
    case 0xbf60c34ef4921000:
      y = 9.9529940556027829412499841055250726640224456787109375e-01;
      break;
    case 0xbf582c1c01d78000:
      y = 9.9660857770005206734964531278819777071475982666015625e-01;
      break;
    case 0x3f5fb29438864810:
      y = 1.00446467107229153725711512379348278045654296875;
      break;
    case 0x3f6d924c938a8000:
      y = 1.0083464751726924912844651771592907607555389404296875;
      break;
    default :
      y = 5.47230882554721187460700093652121722698211669921875e-01;
      y *= R;
      y += 1.1712656335427060749765360014862380921840667724609375;
      y *= R;
      y += 2.03467846365518756357460006256587803363800048828125;
      y *= R;
      y += 2.6509490550833501032457206747494637966156005859375;
      y *= R;
      y += 2.30258509299429903194322832860052585601806640625;
      y *= R;
      y += 1.0;
  }
  
  // Perform output compensation
  y *= ldexp(exp2JBy64[J], M);
  return y;
}
