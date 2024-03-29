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
#include "sinhcosh.h"

#define CONST64BYLN2 92.332482616893656768297660164535045623779296875
#define LN2BY64 0.01083042469624914509729318723429969395510852336883544921875

double rlibm_coshf(float x) {
    float_x fx;
    fx.f = x;
    fx.x &= 0x7FFFFFFF;
    
    // Take care of special cases
    if (fx.x <= 973078527) {
      if (fx.x == 0) return 1.0;
      if (fx.x <= 968164595) return 1.0000000298023223876953125;
      return 1.0000000894069671630859375;
    }
    
    if (fx.x >= 1119016189) {
        if (fx.x > 0x7F800000) return 0.0/0.0;
	if (fx.x == 0x7F800000) return 1.0 / 0.0;
        return 3.40282361850336062550457001444955389952e+38;
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
  
    // Compute sinh component
    double sinhL = 8.33625259341946346636209597136257798410952091217041015625e-03;
    sinhL *= R2;
    sinhL += 1.666666664024638866425931382764247246086597442626953125e-01;
    sinhL *= R2;
    sinhL += 1.0000000000000017763568394002504646778106689453125;
    sinhL *= R;
    
    // Compute cosh component
    double coshL = 4.1668689285602612815129219825394102372229099273681640625e-02;
    coshL *= R2;
    coshL += 4.99999999833178832009394909619004465639591217041015625e-01;
    coshL *= R2;
    coshL += 1.000000000000000444089209850062616169452667236328125;
    
    // Perform output compensation
    return sinhHM * sinhL + coshHM * coshL;
}
