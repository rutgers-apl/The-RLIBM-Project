#include "rlibm.h"

#include <math.h>
#include <mpfr.h>

double getFloat34RNO(mpfr_t mval, mpfr_t r, int sticky) {
  double_x y;
  double max_value = 0x1.ffffff8p+127;
  double max_value_1 = 0x1.ffffffp+127;
  /* Handle case where result is nan */
  if (mpfr_nan_p(mval) != 0) return 0.0f/0.0f;
  /* Handle case where result is zero. If sticky bit is set, return the smallest 34 bit FP number. Else, return 0. */
  if (mpfr_cmp_d(mval, 0.0f) == 0) {
    if (sticky == 0) {
      return 0.0f;
    } else if (sticky < 0) {
      return ldexp(1, -151);
    } else {
      return ldexp(-1, -151);
    }
  }
  if (mpfr_inf_p(mval) && (mpfr_sgn(mval) > 0)) return INFINITY; 
  if (mpfr_inf_p(mval) && (mpfr_sgn(mval) < 0)) return -INFINITY; 
  if (mpfr_cmp_d(mval, max_value_1) > 0) return max_value;
  if (mpfr_cmp_d(mval, -max_value_1) < 0) return -max_value;
  /* Compare mval against the two normal numbers with the normal numbers with the smallest magnitudes to determine if it is denormal. */
  bool mval_is_positive_subnormal = (mpfr_sgn(mval) > 0) && (mpfr_cmp_d(mval, ldexp(1, -126)) < 0);
  bool mval_is_negative_subnormal = (mpfr_sgn(mval) < 0) && (mpfr_cmp_d(mval, ldexp(-1, -126)) > 0);
  if (mval_is_positive_subnormal || mval_is_negative_subnormal) {
    /* mval is a denormal value. */
    /* If the absolute value of mval is less than the second smallest non-zero representable number, just return the number with the smallest representable absolute value (at this point, the number is not zero) */
    if ((mpfr_sgn(mval) > 0) && (mpfr_cmp_d(mval, ldexp(1, -150)) < 0)) return ldexp(1, -151);
    if ((mpfr_sgn(mval) < 0) && (mpfr_cmp_d(mval, ldexp(-1, -150)) > 0)) return ldexp(-1, -151);
    /* Past this point, the most significant bit is not the bit that needs to be updated via the sticky bit. */
    
    mpfr_set_prec(r, 25);
    sticky |= mpfr_set(r, mval, MPFR_RNDZ);
    /* The original mantissa bits will be normalized. */
    y.d = mpfr_get_d(r, MPFR_RNDZ);
    unsigned long exp = y.x >> 52L;
    exp &= 0x7FF;
    exp -= 1023;
    unsigned long shift = -exp - 99;
    if (sticky != 0) y.x |= (1L << shift);
    return y.d;
  }
  mpfr_set_prec(r, 26);
  sticky |= mpfr_set(r, mval, MPFR_RNDZ);
  y.d = mpfr_get_d(r, MPFR_RNDZ);
  /* If there is a sticky bit and the 25th bit after the binary point is 0, it will be updated to 1.
   * If there is a sticky bit and the 25th bit after the binary point is 1, it will remain as 1. */
  if (sticky != 0) y.x |= (1L << 27L);
  return y.d;
}
