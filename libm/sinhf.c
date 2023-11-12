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


#include "rlibm.h"
#include "sinhcosh.h"
#include <math.h>

#define CONST64BYLN2 0x1.71547652b82fep+6
#define LN2BY64 0x1.62e42fefa39efp-7

double rlibm_sinhf(float x) {
  float_x fx;
  fx.f = x;
  unsigned long sign = (fx.x & 0x80000000) == 0 ? 0x0 : 0x8000000000000000;
  fx.x &= 0x7FFFFFFF;

  if (fx.x == 0) return x;
  
  // Take care of special cases
  if (fx.x <= 971544424) {
    double_x dX;
    dX.d = (double)fx.f;
    long exp = (dX.x & 0x7FF0000000000000UL) >> 52UL;
    exp -= 1023L;
    long mantissaCount = exp + 149L;
    if (mantissaCount > 23) mantissaCount = 23;
    mantissaCount += 2L;
    unsigned long shiftAmount = (52L - mantissaCount);
    unsigned long sticky = 1UL << shiftAmount;
    dX.x |= sticky;
    dX.x |= sign;
    return dX.d;
  }
  
  if (fx.x >= 1119016189) {
    if (fx.x > 0x7F800000) return 0.0/0.0;
    if (fx.x == 0x7F800000) {
      if (x > 0.0f) return 1.0 / 0.0;
      else return -1.0 / 0.0;
    }

    if (x > 0.0f) return 0x1.ffffff8p+127;
    else return -0x1.ffffff8p+127;
  }
  
  // Perform range reduction
  double xp = fx.f * CONST64BYLN2;
  int N = (int)xp;
  int N2 = N % 64;
  if (N2 < 0) N2 += 64;
  int N1 = N - N2;
  int I = N1 / 64;
  double R = fx.f - N * LN2BY64;
  double R2 = R * R;
  
  double sinhHigh = sinhKLn2[I];
  double coshHigh = coshKLn2[I];
  double sinhMid = sinhKLn2By64[N2];
  double coshMid = coshKLn2By64[N2];
  
  double sinhHM = sinhHigh * coshMid + coshHigh * sinhMid;
  double coshHM = sinhHigh * sinhMid + coshHigh * coshMid;
  
  // Compute sinh  and coshL component
  double sinhL;
  double coshL;
                           
  if(__builtin_expect(R == 0x1.113e28d466p-7, 0)){
    sinhL = 0x1.113ef7cf95d0cp-7;
    coshL = 0x1.000246c8ff9eep+0;
  }
  else {
    /* sinhL =0x1.ffffffffffffep-1 x^(1) + 0x1.55555554d50dap-3 x^(3) + 0x1.1111b01851046p-7 x^(5) */
    double temp1 = fma(R2,  0x1.1111b01851046p-7,0x1.55555554d50dap-3);
    double temp2 = fma(R2, temp1, 0x1.ffffffffffffep-1);
    sinhL = R * temp2;

    /* coshL =0x1p+0 x^(0) + 0x1.ffffffff997cp-2 x^(2) + 0x1.5555e5da9087ap-5 x^(4) */
    double temp3 = fma (R2, 0x1.5555e5da9087ap-5, 0x1.ffffffff997cp-2);
    coshL = fma(R2, temp3, 0x1p+0);    
  }

  // Perform output compensation
  double_x dX;
  dX.d = sinhHM * coshL + coshHM * sinhL;
  dX.x |= sign;
  return dX.d;
}
