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

#include "rlibm.h"
#include "sinpicospi.h"
#include "pi.h"
#include <math.h>
#include <stdbool.h>

#define PI 0x1.921fb54442d18p+1
#define PI_2 0x1.921fb54442d18p+0

// This is acos.c

double rlibm_acosf(float x) {
  float_x fX;
  fX.f = x;
  fX.x &= 0x7FFFFFFF; 
  if (fX.x > 0x3F800000) {
    return 0.0f/0.0f;
  } 
  if (fX.x < 0x328885A4) {
    return 0x1.921fb58p+0;
  }
  if ((x > 0) && (fX.x == 0x39826222)) {
    return 0x1.920f698p+0;
  }
  double y;
  double R = fX.f;
  R = sqrt(1.0 - R*R)/R;
  bool reciprocal = false;
  if (R > 1.0) {
    R = 1/R;
    reciprocal = true;
  }
  double atan_b = 0.0;
  if (R > 0.001953125) {
    int r;
    double_x dX;
    dX.d = R;
    dX.x -= 1;
    uint32_t value = 0x80 | ((dX.x >> 45) & 0x7f);
    int exponent = dX.x >> 52;
    r = value >> (8 - (exponent - 0x3f6));
    double b = fma(r, 0.00390625, 0.001953125);
    R = (R - b)/fma(R, b, 1.0L);
    atan_b = atan_vals[r];
  }
  double R2 = R*R, temp1, temp2;
  temp1 = fma(R2, 0x1.9e14ca1a8790bp-3, -0x1.55555591c2c0ap-2);
  temp2 = fma(temp1, R2, 0x1.0000000000005p+0);
  y = fma(temp2, R, atan_b);
  if (reciprocal) {
    y = PI_2 - y;
  }
  if (x < 0) y = PI - y;
  return y;
}
