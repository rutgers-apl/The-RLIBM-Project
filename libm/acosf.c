/*
MIT License

Copyright (c) 2022 Sehyeok Park and Santosh Nagarakatte, Rutgers
Architecture and Programming Languages (RAPL) Group

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
#include "pi.h"
#include "rlibm.h"
#define TRUE 1
#define FALSE 0


#define PI 0x1.921fb54442d18p+1
#define PI_2 0x1.921fb54442d18p+0

double rlibm_acosf(float x) {
  float_x fx;
  fx.f = x;
  fx.x &= 0x7FFFFFFF; 
  if (fx.x > 0x3F800000) {
    return 0.0f/0.0f;
  }
  if (fx.x < 0x328885A4) {
    return 0x1.921fb58p+0;
  }
  if (x < 0 && fx.x == 0x39826222) {
    return 0x1.9230018p+0;
  }
  if (fx.x == 0x39826222) {
    return 0x1.920f698p+0;
  }
  double y;
  double R = fx.f;
  R = sqrt(1.0 - R*R)/R;
  int reciprocal = FALSE;
  if (R > 1.0) {
    R = 1/R;
    reciprocal = TRUE;
  }
  double atan_b = 0.0;
  if (R > 0.001953125) {
    int r;
    {
      double_x dx;
      dx.d = R;
      dx.x -= 1;
      uint32_t value = 0x80 | ((dx.x >> 45) & 0x7f);
      int exponent = dx.x >> 52;
      r = value >> (8 - (exponent - 0x3f6));
    }
    double b = r*0.00390625 + 0.001953125;
    R = (R - b)/(1.0L + b*R);
    atan_b = atan_vals[r];
  }
  double R2 = R*R;
  y = 0x1.9e14ca1a8790bp-3;
  y *= R2;
  y += -0x1.55555591c2c0ap-2;
  y *= R2;
  y += 0x1.0000000000005p+0;
  y *= R;
  y += atan_b;
  if (reciprocal) {
    y = PI_2 - y;
  }
  if (x < 0) y = PI - y;
  return y;
}
