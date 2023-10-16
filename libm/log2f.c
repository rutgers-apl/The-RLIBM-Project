/*
  MIT License

  Copyright (c) 2023 Santosh Nagarakatte, Jay Lim, Sehyeok Park, and
  Mridul Aanjaneya, Rutgers Architecture and Programming Languages
  (RAPL) Group

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sub-license, and/or sell
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

/* 
   Implementation of correctly rounded log2(x) in single precision with
   all four rounding modes (round-to-nearest, round-to-zero, round-up,
   round-down) and for all representations from 10-to-32-bits.

   The range reduction step converts a 32-bit input x to a reduced
   input f as follows.  First, the input x is expressed as as x' *
   2^{exp}. Here, x' is in the range [1, 2).

   If the input is subnormal, it is normalized by multiplying by
   2^{23} and subtracting exp by 23.

   Then, x' is range reduced as follows. The reduced input x' (in the form 1.m)  is
   expressed as x' = f + F = F(f/F + 1). Here, F is the value of the
   number whose highest 8-bits including the hidden 1-bit are
   equivalent to x' and the remaining 16 bits are zeros. Hence, f
   represents the difference (x' - F).

   log2(x) = log2(x' * 2^{exp}) = log2(x') + exp
           = log2(F(f/F + )) + exp
           = log2(F) + log2(f/F + 1) + exp

   We create look up tables for all values of log2(F). There are 128
   entries because the leading bit is always 1. Similarly, we create
   lookup tables with 128 entries for 1/F.

   Finally, log2(f/F + 1) is approximated by a polynomial P(f/F)
   generated with the RLIBM method as described in the papers at
   https://people.cs.rutgers.edu/~sn349/rlibm/
*/

#include "rlibm.h"
#include "log.h"
#include <math.h>

double rlibm_log2f(float x) {

  float_x inp = {.f = x};
  uint32_t ux = inp.x;
  uint64_t m = ux & 0x7FFFFF;
  m = m << 29;
  int exp = (ux >> 23) - 127;
  
  if(__builtin_expect(ux < 0x800000 || ux >= 0x7F800000, 0)){

    /* This code for handling subnormals and special cases is from the
       CORE-MATH project:
       https://gitlab.inria.fr/core-math/core-math/-/blob/master/src/binary32/log2/log2f.c
    */    
    if (ux==0||ux==(1u<<31))
      return -__builtin_inff(); // +0.0 || -0.0

    uint32_t inf_or_nan = ((ux>>23)&0xff) == 0xff, nan = inf_or_nan && (ux<<9);

    if (ux>>31 && !nan) return __builtin_nanf("-");

    if (inf_or_nan) return x;

    // subnormal
    int nz = __builtin_clzll(m);
    m <<= nz-11;
    m &= ~0ul>>12;
    exp = exp - (nz - 12);
  }

  /* power of 2 */
  if(__builtin_expect(!m, 0)) return exp;
  
  double_x  xd = {.x = m | 0x3FF0000000000000ULL};
  uint64_t FIndex = m>> 45;
  uint64_t fm = (FIndex) << 45;
  double_x  xf = {.x = fm |0x3FF0000000000000ULL};
  double f = xd.d - xf.d;  
  
  f *= __log_oneByF[FIndex];

  double coeffs[] = {
    0x1.71547652bcde3p+0,
    -0x1.7154769679dd8p-1,
    0x1.ec7198ec61291p-2,
    -0x1.72033bee9c2d6p-2,
    0x1.4f082e01903edp-2
  };
  
  double xsquare = f*f;

  double temp1 = coeffs[3] + f * coeffs[4];
  double temp2 = coeffs[1] + f * coeffs[2];
  double temp3 = temp2 + xsquare * temp1;
  double temp4 = xsquare * temp3;
  double temp5 = __log2_lut[FIndex] + exp;
  double y = temp4 + f * coeffs[0] + temp5;
  return y;
}
