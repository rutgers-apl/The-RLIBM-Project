/*

MIT License

Copyright (c) 2023 Santosh Nagarakatte, Jay Lim, Sehyeok Park, and
Mridul Aanjaneya, Rutgers Architecture and Programming Languages
(RAPL) Group

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
#include "exp2.h"

double rlibm_expf(float x) { 
  double_x dY;
  float_x fx;
  fx.f = x;

  // Take care of special cases
  if ((fx.x & 0x7FFFFFFF) == 0) return 1.0;

  if (fx.x <= 0x33ffffff) {
    if (fx.x <= 0x337fffff) return 0x1.0000008p+0;
    return 0x1.0000018p+0;
  }

  if (0x42b17218 <= fx.x && fx.x <= 0xb3800000) {
    if (fx.x < 0x80000000) {
      if (fx.x < 0x7F800000) return 0x1.ffffff8p+127;
      if (fx.x == 0x7F800000) return 1.0 / 0.0;
      return 0.0/0.0;
    }
    
    // negative counterpart
    if (fx.x <= 0xb3000000) return 0x1.ffffff8p-1; 
    
    return 0x1.fffffe8p-1; 
  }

  if (fx.x >= 0xc2cff1b5) {
    if (fx.x == 0xFF800000) return 0.0;
    if (fx.x < 0xFF800000) {
      dY.x = 0x3680000000000000;
      return dY.d;
    }
    return 0.0/0.0;
  }
  
  // Perform range reduction
  double xp = x * 0x1.71547652b82fep+6;
  int N = (int)xp;
  int N2 = N % 64;
  if (N2 < 0) N2 += 64;
  int N1 = N - N2;
  
  int M = N1 / 64;
  int J = N2;
  double R = x - N * 0x1.62e42fefa39efp-7;
  
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
  dY.d = exp2JBy64[J];
  dY.x += ((uint64_t)M << 52);
  y *= dY.d;
  return y;
}
