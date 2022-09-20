#include <math.h>
#include "rlibm.h"
#include "exp2.h"

double rlibm_exp(float x) {
  float_x fx;
  fx.f = x;
  
  // Take care of special cases
  if ((fx.x & 0x7FFFFFFF) == 0) return 1.0;

  if (fx.x <= 872415231) {
    if (fx.x <= 864026623) return 1.0000000298023223876953125;
    return 1.0000000894069671630859375;
  }

  if (1118925336 <= fx.x && fx.x <= 3011510272) {
    if (fx.x < 0x80000000) {
      if (fx.x < 0x7F800000) return 3.40282361850336062550457001444955389952e+38;
      if (fx.x == 0x7F800000) return 1.0 / 0.0;
      return 0.0/0.0;
    }
    
    // negative counterpart
    if (fx.x <= 3003121664) return 0.99999998509883880615234375;
    
    return 0.99999995529651641845703125;
  }

  if (fx.x >= 3268407733) {
    if (fx.x == 0xFF800000) return 0.0;
    if (fx.x < 0xFF800000) return ldexp(1.0, -151);
    return 0.0/0.0;
  }
  
  // Perform range reduction
  double xp = x * 92.332482616893656768297660164535045623779296875;
  int N = (int)xp;
  int N2 = N % 64;
  if (N2 < 0) N2 += 64;
  int N1 = N - N2;
  
  int M = N1 / 64;
  int J = N2;
  double R = x - N *
  0.01083042469624914509729318723429969395510852336883544921875;
  
  double_x dX;
  dX.d = R;
  double y;
 
  if (R < -0x1.9e76dcp-24){
    if(R == -0x1.f925ff514p-9){
      y = 0x1.fe07d2e08aadp-1 ;
    } else if(R == -0x1.e08d3ep-9){
      y = 0x1.fe20540000001p-1;
    } else if(R == -0x1.c57f87ab4d8p-9){
      y = 0x1.fe3b49143b1b9p-1;
    } else{
      double coeffs[]  = {
			  0x1.0000000000004p+0,
			  0x1.fffffffffa598p-1,
			  0x1.ffffffaabbeacp-2,
			  0x1.5554722a9f8d2p-3,
			  0x1.53a5b76499a6cp-5
      };

      double xsquare = R * R;
      double temp1 = fma(R, coeffs[1], coeffs[0]);
      double temp2 = fma(R, coeffs[3], coeffs[2]);
      double temp3 = fma(xsquare, coeffs[4], temp2);
      y = fma(xsquare, temp3, temp1);
    }
  } else {

    
    if(R == 0x1.113e28d466p-7){
      y = 0x1.0224c4b89eca9p+0;
    } else if(R == 0x1.5d1488105c611p-7){
      y = 0x1.02bde475d8a85p+0;
    } else if(R == 0x1.5f25f7dc62ap-7){
      y = 0x1.02c212b5454bbp+0;
    }else{

      double coeffs[]  = {
			  0x1.ffffffffffffep-1,
			  0x1.fffffffffac51p-1,
			  0x1.0000001e24d1bp-1,
			  0x1.5554a0a1c255bp-3,
			  0x1.56d939dd19dbfp-5
      };
      double xsquare = R * R;
      double temp1 = fma(R, coeffs[1], coeffs[0]);
      double temp2 = fma(R, coeffs[3], coeffs[2]);
      double temp3 = fma(xsquare, coeffs[4], temp2);
      y = fma(xsquare, temp3, temp1);
    }
  }
  
  // Perform output compensation
  y *= ldexp(exp2JBy64[J], M);
  return y;
}
