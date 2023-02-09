/*
MIT License
Copyright (c) 2023 Sehyeok Park, Mridul Aanjaneya, and Santosh
Nagarakatte, Rutgers Architecture and Programming Languages (RAPL)
Group
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
#include "sinpicospi.h"
#include "pi.h"

#define PI 0x1.921fb54442d18p+1
#define PI_OVER_2 0x1.921fb54442d18p+0

double rlibm_acosf(float x) {
  float_x fX;
  fX.f = x;
  int exp = (fX.x>>23)&0xff;
  if (exp >= 127) {
    if (fX.x == 0x3f800000) return 0.0f;
    if (fX.x == 0xbf800000) return 0x1.921fb58p+1;
    return 0.0f/0.0f;
  }
  if (exp < 126) {
    if (exp < 115) {
      if (fX.x == 0x328885a3) return 0x1.921fb58p+0; 
      return PI_OVER_2 - x;
    }
    if (fX.x == 0x39826222) return 0x1.920f698p+0;
    double R = x;    
    double R2 = R*R, R4 = R2*R2, R8 = R4*R4;
    double temp1 = fma(0x1.5555555697c1bp-3, R2, 0x1.fffffffffff3p-1);
    double temp2 = fma(0x1.6db7bbb6ded7ap-5, R2, 0x1.3333318732c22p-4);
    double temp3 = fma(0x1.736ea9a422144p-6, R2, 0x1.f189f10aca3c4p-6);
    double temp4 = fma(0x1.6e3128df7327fp-5, R2, 0x1.bc662b2b20facp-7);
    double temp5 = fma(0x1.087096701e543p-1, R2, -0x1.2dc923fbefeafp-3);
    double temp6 = fma(0x1.6c756583d7291p-1, R2, -0x1.cd0193b294652p-1);
    temp1 = fma(temp2, R4, temp1);
    temp3 = fma(temp4, R4, temp3);
    temp5 = fma(temp6, R4, temp5);
    temp3 = fma(temp5, R8, temp3);
    temp1 = x*fma(temp3, R8, temp1);
    return PI_OVER_2 - temp1;
  } else {
    double R = 1.0f - fabs(x);
    double s = __builtin_sqrt(R);
    double R2 = R*R, R4 = R2*R2;
    double temp1 = fma(0x1.e2b7dde5f9c6ap-4, R, 0x1.6a09e667f3335p+0);
    double temp2 = fma(0x1.029a9d0b9e527p-7, R, 0x1.b27241ebe50e6p-6);
    double temp3 = fma(0x1.0898e2162443cp-10, R, 0x1.5faf273ebb161p-9);
    double temp4 = fma(0x1.374a8455e9bc5p-11, R, 0x1.15ee663e750d9p-12);
    double temp5 = fma(0x1.a39a22691ed9cp-10, R, -0x1.035c9d5dd6617p-10);
    double temp6 = fma(0x1.0526ef52072a1p-11, R, -0x1.5ad259dd167bfp-10);
    temp1 = fma(temp2, R2, temp1);
    temp3 = fma(temp4, R2, temp3);
    temp5 = fma(temp6, R2, temp5);
    temp3 = fma(temp5, R4, temp3);
    double y = s*fma(temp3, R4, temp1);
    return (fX.x>>31) ? PI - y : y;
  }
}
