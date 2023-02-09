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

double rlibm_atanf(float x) {
  double y;
  float_x fX;
  fX.f = x;
  fX.x &= 0x7FFFFFFF;
  double R, R2;
  if (fX.x < 0x39B89BA3) {
    return 0x1.ffffffffffffep-1*x;
  } else if (fX.x <= 0x4C700517) {
    if (fX.x == 0x3AAC434B) {
      return __builtin_copysign(0x1.5886888p-10, x);
    }
    if (fX.x == 0x3FEEFCFB) {
      return __builtin_copysign(0x1.143ec48p+0, x);
    }
    R = x;
    int gt = ((fX.x>>23)>=127);
    if (gt) R=1/R;
    int idx = (int)(fabs(R*256.0f));
    double atan_b = atan_vals_256[idx];
    if (idx) {
      double_x b;
      b.d = R;
      uint64_t exp = ((b.x>>52)&0x7ff);
      b.x = ((b.x>>(1067-exp))<<(1067-exp));
      R = (R-b.d)/(R*b.d + 1.0f);
    }
    R2 = R*R;
    double temp1 = fma(0x1.99af64261375p-3, R2, -0x1.55555550d9454p-2);
    double temp2 = fma(0x1.0000000000002p+0, R, __builtin_copysign(atan_b,x));
    y = fma(temp1, R2*R, temp2);

    if (gt) {
      return __builtin_copysign(PI_OVER_2,x) - y;
    } else {
      return y;
    }
  } else if (fX.x <= 0x7F800000) {
    return __builtin_copysign(0x1.921fb5fffffffp+0, x);
  } else {
    return 0.0f/0.0f;
  }
}
