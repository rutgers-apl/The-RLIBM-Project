/*

MIT License

Copyright (c) 2022 Santosh Nagarakatte, Jay Lim, Sehyeok Park, and
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


#include "rlibm.h"
#include "log.h"
#include <math.h>

#define LN2HIGH 0.69314718055994528622676398299518041312694549560546875

double rlibm_logf(float x) {
  float_x fix, fit;
  fix.f = x;
  int m = 0;
  
  if (x == 1.0) return 0.0;
  
  if (fix.x < 0x800000 || fix.x >= 0x7F800000) {
    if ((fix.x & 0x7FFFFFFF) == 0) { // log(+/-0) = -infty
      fix.x = 0xFF800000;
      return fix.f;
    }
    
    if (fix.x > 0x7FFFFFFF) return (x - x) / 0; // Log(-val) = NaN
    if (fix.x >= 0x7F800000) return x + x;
    fix.f *= 8.388608e+06;
    m -= 23;
  }
  
  m += fix.x >> 23;
  m -= 127;
  fix.x &= 0x007FFFFF;
  fix.x |= 0x3F800000;
  
  fit.x = fix.x & 0x007F0000;
  int FIndex = fit.x >> 16;
  fit.x |= 0x3F800000;
  
  double f = fix.f - fit.f;
  f *= __log_oneByF[FIndex];
  
  double y = 0.0;
  
  if (f < 0x1.5a8f8d28ac42fp-9) {

    if(f == 0x1.67f6db6db6db7p-10) {
      y = 0x1.67b7a57460101p-10;
    } else if(f == 0x1.d6f7e432f7e44p-10){
      y = 0x1.d68bb6f7c2101p-10;
    } else if(f == 0x1.ebdfbf7efdfbfp-10){
      y = 0x1.eb69c2c46f5p-10;
    } else if(f == 0x1.fbd2361d2361ep-10){
      y = 0x1.fb54746536101p-10;
    } else if(f == 0x1.57497b425ed09p-9){
      y = 0x1.56d6991a29e81p-9;
    } else{

      double coeffs[4]={
			0x1.fffffffffffd1p-1,
			-0x1.ffffffa75f811p-2,
			0x1.55539dcca455cp-2,
			-0x1.fb54dd294958dp-3		
      };
      /* polynomial is y=0x1.fffffffffffd1p-1 x^(1) +
	 -0x1.ffffffa75f811p-2 x^(2) + 0x1.55539dcca455cp-2 x^(3) +
	 -0x1.fb54dd294958dp-3 x^(4) */

      double temp1 = fma(f, coeffs[2], coeffs[1]);
      double xsquare = f*f;
      double temp2 = fma(xsquare, coeffs[3], temp1);
      double temp3 = xsquare * temp2;
      y = fma(f, coeffs[0], temp3);
    }
  } else {

    if(f == 0x1.78d3dcb08d3ddp-9){
      y = 0x1.78496eb2b2bc1p-9;      
    } else if(f == 0x1.e8a1fd1b7af01p-9){
      y = 0x1.e7b9668c47041p-9;      
    } else if(f == 0x1.3155555555555p-8){
      y = 0x1.309fcf6432ca1p-8;      
    } else if(f == 0x1.740a7ac29eb0ap-8){
      y = 0x1.72fd098c2a761p-8;      
    } else if(f == 0x1.a486d6f63aa04p-8){
      y = 0x1.a32ee8debeb9bp-8;
    } else if(f == 0x1.f9fcp-8){
      y = 0x1.f80a850000001p-8;
    } else{
      double coeffs[4]={
			0x1.ffffffff5adf6p-1,
			-0x1.fffffb33f4724p-2,
			0x1.554ed3225f9f3p-2,
			-0x1.f869f77cf2e1cp-3
      };
      /* Polynomial: y=0x1.ffffffff4dc8bp-1 x^(1) +
	 -0x1.fffffaeece64cp-2 x^(2) + 0x1.554e9b388c19ep-2 x^(3) +
	 -0x1.f84d9fdd570dp-3 x^(4) */

      double temp1 = fma(f, coeffs[2], coeffs[1]);
      double xsquare = f*f;
      double temp2 = fma(xsquare, coeffs[3], temp1);
      double temp3 = xsquare * temp2;
      y = fma(f, coeffs[0], temp3);
    }    
  }
  
  y += __ln_lutHIGH[FIndex];
  y += m * LN2HIGH;
  
  return y;
}
