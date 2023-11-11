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

/* 
   Implementation of correctly rounded log10f(x) in single precision with
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

   log10f(x) = log10f(x' * 2^{exp}) = log10f(x') + exp * log10 2
           = log10f(F(f/F + 1)) + exp * log10 2
           = log10f(F) + log10f(f/F + 1) + exp

   We create look up tables for all values of log10f(F). There are 128
   entries because the leading bit is always 1. Similarly, we create
   lookup tables with 128 entries for 1/F.

   Finally, log10f(f/F + 1) is approximated by a polynomial f* P(f/F) 
   generated with the RLIBM method as described in the papers at
   https://people.cs.rutgers.edu/~sn349/rlibm/
*/



#include "rlibm.h"
#include "log.h"
#include <math.h>
#include <assert.h>

#define LOG102HIGH 0x1.34413509f79fep-2
#define LOG102LOW  0x1.e623e2566b02ep-55


double rlibm_log10f(float x) {

  //  double old_rr = old_range_reduction(x);
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

  double_x  xd = {.x = m | 0x3FF0000000000000ULL};
  uint64_t FIndex = m>> 45;
  uint64_t fm = (FIndex) << 45;
  double_x  xf = {.x = fm |0x3FF0000000000000ULL};
  double f = xd.d - xf.d;  

  f *= __log_oneByF[FIndex];

  double coeffs[6] =
    {
      0x1.bcb7b1526eba1p-2,
      -0x1.bcb7b152ff289p-3,
      0x1.287a6320ef253p-3,
      -0x1.bc614c3bdbca8p-4,
      0x1.23bf152c3f48ap-4,
      0x1.cb2f0b48e857fp-1,
    };
  double y = exp * LOG102LOW + _rlibm_lut_log10F[FIndex] + exp * LOG102HIGH;

  if(__builtin_expect(f == 0x1.bde34a2b10bf6p-9,  0)) {
    return 0x1.82a2d0cd6ee1p-10  + y;
  }
  
  if(__builtin_expect(f == 0x1.8ff099fc267fp-9,  0)) {
    return  0x1.5adab7bb93889p-10 + y;
  }


  if(__builtin_expect(f ==0x1.8abe0f83e0f84p-9,  0)) {
    return 0x1.565a88a764385p-10  + y;
  }

  if(__builtin_expect(f == 0x1.8533a8c0dc69ap-10 ,  0)) {
    return 0x1.51ce40d6f4f67p-11  + y;
  }

  if(__builtin_expect(f == 0x1.2bec04fec04ffp-8,  0)) {
    return  0x1.03ea0f7b7496dp-9  + exp * LOG102LOW + _rlibm_lut_log10F[FIndex] + exp * LOG102HIGH;
  }
  
  if(__builtin_expect(f == 0x1.2af84a062b2e5p-8,  0)) {
    return 0x1.03176656347c4p-9 + exp * LOG102LOW + _rlibm_lut_log10F[FIndex] + exp * LOG102HIGH;
  }


  if(__builtin_expect(f == 0x1.e0b56ad5ab56bp-9,  0)) {
    return 0x1.a0c640ef6a074p-10 + y;
  }

  if(__builtin_expect(f == 0x1.1fddb0d3224f3p-9,  0)) {
    return 0x1.f386956531508p-11 + exp * LOG102LOW + _rlibm_lut_log10F[FIndex] + exp * LOG102HIGH;
  }

  double xsquare = f * f;
  double xcube = f * xsquare;
  double temp1 = fma(coeffs[1] ,f , coeffs[0]);
  double temp2 = fma(coeffs[2], xsquare, temp1);
  double temp3 = fma(coeffs[4], f , coeffs[3]);
  double temp4 = fma(coeffs[5], xsquare, temp3);
  double temp5 = fma(temp4, xcube ,  temp2);
  double temp6 =  f * temp5;
  return temp6 + exp * LOG102LOW + _rlibm_lut_log10F[FIndex] + exp * LOG102HIGH;  
}
