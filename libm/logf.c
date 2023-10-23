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
   Implementation of correctly rounded logf(x) in single precision with
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

   logf(x) = logf(x' * 2^{exp}) = logf(x') + exp
           = logf(F(f/F + 1)) + exp
           = logf(F) + logf(f/F + 1) + exp

   We create look up tables for all values of logf(F). There are 128
   entries because the leading bit is always 1. Similarly, we create
   lookup tables with 128 entries for 1/F.

   Finally, logf(f/F + 1) is approximated by a polynomial f* P(f/F) 
   generated with the RLIBM method as described in the papers at
   https://people.cs.rutgers.edu/~sn349/rlibm/
*/



#include "rlibm.h"
#include "log.h"
#include <math.h>
#include <assert.h>

#define LN2HIGH  0x1.62e42fefa39efp-1

double rlibm_logf(float x) {

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


  double coeffs[5] = {
    0x1.0000000007a9cp+0,
    -0x1.0000004ef6ca5p-1,
    0x1.55565c6f85cc2p-2,
    -0x1.00af53c746b32p-2,    
    0x1.e767c61247c1ep-3
  };

  if(__builtin_expect(f == 0x1.7096969696969p-11, 0)) {
      return 0x1.707567c76c101p-11 + exp * LN2HIGH + rlibm_lnF[FIndex]; 
  }

  if(__builtin_expect(f == 0x1.67f6db6db6db7p-10,  0)) {
      return 0x1.67b7a57462001p-10 + exp * LN2HIGH + rlibm_lnF[FIndex]; 
  }

  if(__builtin_expect(f == 0x1.d6f7e432f7e44p-10, 0)) {
    return 0x1.d68bb6f7c2101p-10 +  exp * LN2HIGH + rlibm_lnF[FIndex]; 
  }

  if(__builtin_expect(f == 0x1.f96p-10, 0)) {
    return 0x1.f8e3730000001p-10 +  exp * LN2HIGH + rlibm_lnF[FIndex]; 
  }

  if(__builtin_expect(f == 0x1.fbd2361d2361ep-10, 0)) {
    return 0x1.fb54746536101p-10 +  exp * LN2HIGH + rlibm_lnF[FIndex]; 
  }

  
  if(__builtin_expect(f == 0x1.23624dd2f1aap-9, 0)) {
    return 0x1.230f831236001p-9  + exp * LN2HIGH + rlibm_lnF[FIndex]; 
  }

  if(__builtin_expect(f == 0x1.57497b425ed09p-9, 0)) {
    return 0x1.56d6991a2a001p-9 +  exp * LN2HIGH + rlibm_lnF[FIndex]; 
  }

  if(__builtin_expect(f == 0x1.bab0df6b0df6bp-9, 0)) {
    return 0x1.b9f1e20cc6801p-9  +  exp * LN2HIGH + rlibm_lnF[FIndex]; 
  }

  if(__builtin_expect(f == 0x1.e8a1fd1b7af01p-9, 0)) {
    return 0x1.e7b9668c44001p-9 + exp * LN2HIGH + rlibm_lnF[FIndex]; 
  }
  if(__builtin_expect(f == 0x1.05c77b03531dfp-8, 0)) {
    return 0x1.0541fbde4def1p-8 + exp * LN2HIGH + rlibm_lnF[FIndex]; 
  }
  
  if(__builtin_expect(f == 0x1.3155555555555p-8, 0)) {
    return 0x1.309fcf6433001p-8 + exp * LN2HIGH + rlibm_lnF[FIndex]; 
  }

  
  if(__builtin_expect(f == 0x1.9c0a2c145828bp-8, 0)) {
    return 0x1.9abff5d8ca001p-8 +  exp * LN2HIGH + rlibm_lnF[FIndex];
  }

  if(__builtin_expect(f == 0x1.a33ep-8, 0)) {
    return 0x1.a1e82b0000001p-8 +  exp * LN2HIGH + rlibm_lnF[FIndex];
  }

  
  if(__builtin_expect(f == 0x1.c15d1745d1746p-8, 0)) {
    return 0x1.bfd47e987b801p-8  + exp * LN2HIGH + rlibm_lnF[FIndex];
  }


  double temp1 = f * coeffs[4];
  double temp2 = f * (coeffs[3] + temp1);
  double temp3 = f * (coeffs[2] + temp2);
  double temp4 = f * (coeffs[1] + temp3);
  double temp5 = f * (coeffs[0] + temp4);

  double y = exp * LN2HIGH + rlibm_lnF[FIndex];
  return temp5 +  y;
  
}
