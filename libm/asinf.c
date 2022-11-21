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

#define PI_2 0x1.921fb54442d18p+0

double rlibm_asinf(float x) {
  double_x dY;
  float_x fX;
  fX.f = x;
  fX.x &= 0x7FFFFFFF;
  uint64_t sgn = (x < 0) ? 0x8000000000000000 : 0; 
  if (fX.x > 0x3F800000) {
    return 0.0f/0.0f;
  }
  double R, R2;
  if (fX.x < 0x39E89768) {
    dY.d = 0x1.0000008000007p+0*fX.f;
    dY.x ^= sgn;
    return dY.d;
  } else if (fX.x <= 0x3B000000) {
    R = fX.f;
    R2 = R*R;
    double temp = fma(R2, 0x1.40076be912c25p-4, 0x1.5555534e4e023p-3);
    dY.d = R*fma(temp, R2, 0x1.000000000003fp+0);
    dY.x ^= sgn;
    return dY.d;
  } else {
    if (fX.x == 0x3B5637DC) {
      dY.d = 0x1.ac6fe98p-9;
      dY.x ^= sgn;
      return dY.d;
    }
    if (fX.x == 0x3B7A8D0D) {
      dY.d = 0x1.f51a698p-9;
      dY.x ^= sgn;
      return dY.d;
    }
    if (fX.x == 0x3F083A1A) {
      dY.d = 0x1.1f4b648p-1;
      dY.x ^= sgn;
      return dY.d;
    }
    R = fX.f;
    R /= sqrt(1.0 - R*R);
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
    dY.d = fma(temp2, R, atan_b);
    if (reciprocal) {
      dY.d = PI_2 - dY.d;
    }
  }
  dY.x ^= sgn;
  return dY.d;
}
