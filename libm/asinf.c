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

#define PI_OVER_2 0x1.921fb54442d18p+0

double rlibm_asinf(float x) {
  float_x fX;
  fX.f = x;
  int exp = ((fX.x>>23)&0xff);
  if (exp < 126) {
    if (exp < 115) return 0x1.0000008000007p+0*x;
    double R = x;
    double R2 = R*R, R4 = R2*R2, R8 = R4*R4;
    double temp1 = fma(0x1.55555555493abp-3, R2, 0x1.0000000000001p+0);
    double temp2 = fma(0x1.6db6c8070fa5ep-5, R2, 0x1.3333334def48dp-4);
    double temp3 = fma(0x1.6e15ba5571e71p-6, R2, 0x1.f1ccfa53189cdp-6);
    double temp4 = fma(0x1.7d6d85305fe3bp-7, R2, 0x1.21ad99171ef4p-6);
    double temp5 = fma(-0x1.f35cff48572ccp-7, R2, 0x1.61b07ca56ea33p-6);
    temp1 = fma(temp2, R4, temp1); 
    temp3 = fma(temp4, R4, temp3); 
    temp5 = fma(0x1.3ca9d85592ce8p-5, R4, temp5); 
    temp3 = fma(temp5, R8, temp3);
    return x*fma(temp3, R8, temp1);
  } else if (exp == 126) {
    double R = 1.0f - fabs(x);
    double s = __builtin_sqrt(R);
    double R2 = R*R, R4 = R2*R2;
    double temp1 = fma(0x1.e2b7ddf9799c8p-4, R, 0x1.6a09e667efab3p+0);
    double temp2 = fma(0x1.029bef197bd55p-7, R, 0x1.b27237cd13384p-6);
    double temp3 = fma(0x1.0ac5b654d7e7p-10, R, 0x1.5f7d64a1619e8p-9);
    double temp4 = fma(0x1.50cf0a55e80ffp-11, R, 0x1.f231c86dd72f9p-13);
    double temp5 = fma(0x1.8f4b1c200396ap-10, R, -0x1.08efc3da631b9p-10);
    double temp6 = fma(0x1.bf695e7e14cb6p-12, R, -0x1.36b03676b2a27p-10);
    temp1 = fma(temp2, R2, temp1);
    temp3 = fma(temp4, R2, temp3);
    temp5 = fma(temp6, R2, temp5);
    temp3 = fma(temp5, R4, temp3);
    temp1 = fma(temp3, R4, temp1);
    double y = fma(-s, temp1, PI_OVER_2);
    return __builtin_copysign(y, x);
  } else {
    if ((fX.x&0x7fffffff) == 0x3f800000) return __builtin_copysign(0x1.921fb58p+0,x);
    return 0.0/0.0f;
  }
}
