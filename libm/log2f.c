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
#include "log2.h"
#include <math.h>

double rlibm_log2f(float x) {
  float_x fix, fit, spec;
  fix.f = x;
  int m = 0;

  if (fix.x < 0x800000 || fix.x >= 0x7F800000) {
    if ((fix.x & 0x7FFFFFFF) == 0) { // log(+/-0) = -infty
      fix.x = 0xFF800000;
      return fix.f;
    }
    
    if (fix.x > 0x7FFFFFFF) { // Log(-val) = NaN
      return (x - x) / 0;
        
    }
    
    if (fix.x >= 0x7F800000) {
      return x + x;
    }

    // Special case when we have denormal input and exact result
    int exp;
    spec.f = frexpf(fix.f, &exp);
    if (spec.x == 0x3f000000) return (double)(exp - 1);

    fix.f *= 8.388608e+06;
    m -= 23;
  }
  
  m += fix.x >> 23;
  m -= 127;
  fix.x &= 0x007FFFFF;
  
  if (fix.x == 0) {
    return (double)m;
  }
  
  fix.x |= 0x3F800000;
  
  fit.x = fix.x & 0x007F0000;
  int FIndex = fit.x >> 16;
  fit.x |= 0x3F800000;
  
  double f = fix.f - fit.f;
  f *= __log_oneByF[FIndex];

  double y = 0.0;


  // FMA coefficients
  double coeffs[] = {
		     1.4426950408932206482148785653407685458660125732421875000000000000000000e+00,
		     -7.2134752833238291458428648184053599834442138671875000000000000000000000e-01,
		     4.8090208952661678276641055163054261356592178344726562500000000000000000e-01,
		     -3.6134046214924120388189976438297890126705169677734375000000000000000000e-01,
		     3.2717964062254517587646773790766019374132156372070312500000000000000000e-01		     
  };
  
  double xsquare = f*f;
  double temp1 = fma(f, coeffs[4], coeffs[3]);
  double temp2 = fma(f, coeffs[2], coeffs[1]);
  double temp3 = fma(xsquare, temp1, temp2);
  double temp4 = xsquare * temp3;
  y = fma(f, coeffs[0], temp4);    
  y += __log2_lut[FIndex];
  y += m;
  
  return y;
}
