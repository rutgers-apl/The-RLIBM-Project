#pragma once
#include "mpfr.h"
#include "math.h"

// MAXVAL = 0111 1111 0111 1111 1111 1111 1111 1111 11
#define MAXVAL 3.40282361850336062550457001444955389952e+38
// MAXm1VAL = 0111 1111 0111 1111 1111 1111 1111 1111 10
#define MAXm1VAL 3.40282356779733661637539395458142568448e+38
#define OVER34VAL 3.40282366920938463463374607431768211456e+38
#define MINVAL 3.503246160812042677309323958224790328200654854691289429392670709724477706714651503716595470905303955078125e-46

typedef union {
    float f;
    unsigned x;
} floatX;

typedef union {
    double d;
    unsigned long long int x;
} doubleX;

double FromMPFRToFloat34Ro(mpfr_t mval, int sticky) {
    int cmp, status;
    double retVal, ldexpVal;
    
    // Check for really special cases
    if (mpfr_nan_p(mval)) return 0.0/0.0;
    if (mpfr_zero_p(mval) && sticky == 0) return 0.0;
    if (mpfr_inf_p(mval)) return mpfr_get_d(mval, MPFR_RNDZ);

    // Take care of positive largest value:
    cmp = mpfr_cmp_d(mval, MAXm1VAL);
    if (cmp > 0) return MAXVAL;
    if (cmp == 0 && sticky == 1) return MAXVAL;
    
    // Take care of negative largest value:
    cmp = mpfr_cmp_d(mval, -MAXm1VAL);
    if (cmp < 0) return -MAXVAL;
    if (cmp == 0 && sticky == 1) return -MAXVAL;
    
    // Check for denormalized values
    if (mpfr_cmp_d(mval, ldexp(-1.0, -126)) > 0 &&
        mpfr_cmp_d(mval, ldexp(1.0, -126)) < 0) {
    
        if (mpfr_cmp_d(mval, 0.0) > 0) {
            // Value is positive
            
            // Take care of denormalized values:
            // 0000 0000 0000 0000 0000 0000 0000 0000 00 = 0
            // 0000 0000 0000 0000 0000 0000 0000 0000 01 = 2^-151          1
            
            // From here on out, we can use a pattern.
            // 0000 0000 0000 0000 0000 0000 0000 0000 10 = 2^-150          2
            // 0000 0000 0000 0000 0000 0000 0000 0000 11 = 1.5 * 2^-150    2
            // 0000 0000 0000 0000 0000 0000 0000 0001 00 = 2^-149          3
            // 0000 0000 0000 0000 0000 0000 0000 0010 00 = 2^-148          4
            // 0000 0000 0000 0000 0000 0000 0000 0100 00 = 2^-147          5
            // 0000 0000 0000 0000 0000 0000 0000 1000 00 = 2^-146          6
            // 0000 0000 0000 0000 0000 0000 0001 0000 00 = 2^-145          7
            // 0000 0000 0000 0000 0000 0000 0010 0000 00 = 2^-144          8
            // 0000 0000 0000 0000 0000 0000 0100 0000 00 = 2^-143          9
            // 0000 0000 0000 0000 0000 0000 1000 0000 00 = 2^-142          10
            // 0000 0000 0000 0000 0000 0001 0000 0000 00 = 2^-141          11
            // 0000 0000 0000 0000 0000 0010 0000 0000 00 = 2^-140          12
            // 0000 0000 0000 0000 0000 0100 0000 0000 00 = 2^-139          13
            // 0000 0000 0000 0000 0000 1000 0000 0000 00 = 2^-138          14
            // 0000 0000 0000 0000 0001 0000 0000 0000 00 = 2^-137          15
            // 0000 0000 0000 0000 0010 0000 0000 0000 00 = 2^-136          16
            // 0000 0000 0000 0000 0100 0000 0000 0000 00 = 2^-135          17
            // 0000 0000 0000 0000 1000 0000 0000 0000 00 = 2^-134          18
            // 0000 0000 0000 0001 0000 0000 0000 0000 00 = 2^-133          19
            // 0000 0000 0000 0010 0000 0000 0000 0000 00 = 2^-132          20
            // 0000 0000 0000 0100 0000 0000 0000 0000 00 = 2^-131          21
            // 0000 0000 0000 1000 0000 0000 0000 0000 00 = 2^-130          22
            // 0000 0000 0001 0000 0000 0000 0000 0000 00 = 2^-129          23
            // 0000 0000 0010 0000 0000 0000 0000 0000 00 = 2^-128          24
            // 0000 0000 0100 0000 0000 0000 0000 0000 00 = 2^-127          25
            // # precision = (152 + exp)
            // First normalized value:
            // 0000 0000 1000 0000 0000 0000 0000 0000 00 = 2^-126          26
            
            
            // 0000 0000 0000 0000 0000 0000 0000 0000 01 = 2^-151
            // 0000 0000 0000 0000 0000 0000 0000 0000 10 = 2^-150
            // if _mval < 2^-150, then round to 2^-151.
            ldexpVal = ldexp(1.0, -150);
            retVal = ldexp(1.0, -151);
            cmp = mpfr_cmp_d(mval, ldexpVal);
            if (cmp < 0) return retVal;
            
            // 0000 0000 0000 0000 0000 0000 0000 0000 10 = 2^-150
            // if _mval == 2^-150 and sticky == 0, then it's exact: 2^-150
            retVal = ldexpVal;
            if (cmp == 0 && sticky == 0) return retVal;
            
            // 0000 0000 0000 0000 0000 0000 0000 0000 10 = 2^-150
            // 0000 0000 0000 0000 0000 0000 0000 0000 11 = 1.5 * 2^-150
            // 0000 0000 0000 0000 0000 0000 0000 0001 00 = 2^-149
            // if _mval < 2^-149, then round to 1.5 * 2^-150
            ldexpVal = ldexp(1.0, -149);
            retVal = ldexp(1.5, -150);
            cmp = mpfr_cmp_d(mval, ldexpVal);
            if (cmp < 0) return retVal;
            
            // 0000 0000 0000 0000 0000 0000 0000 0001 00 = 2^-149
            // if _mval == 2^-149 and sticky == 0, then it's exact: 2^-149
            retVal = ldexpVal;
            if (cmp == 0 && sticky == 0) return retVal;
            
            // From here on out, we have at least 3 precision bits, so we can work
            // with the normal logic
            
            // 1. Find out how many precision bits we will have
            long exp;
            double fr = mpfr_get_d_2exp(&exp, mval, MPFR_RNDZ);
            long origExp = exp;
            fr *= 2;
            exp--;
            
            // 2. Round _mval to a value with (numPrec - 1) precision bits using
            // round towards zero. Keep track of sticky bit
            long numPrec = 152 + exp;
            mpfr_t r;
            mpfr_init2(r, numPrec - 1);
            status = mpfr_set(r, mval, MPFR_RNDZ);
            if (status != 0) sticky |= 0x1;
            
            // 3. Convert "r" to double. It should be exactly representable.
            doubleX dx;
            dx.d = mpfr_get_d(r, MPFR_RNDZ);
            
            // 4. Add "sticky bit" to the numPrec-th bit
            // Let's say we have 13 precision bits.
            // seee eeee eeee xxxx xxxx xxxx ____ ____ ____ ____ ____ ____ ____ ____ ____ ____
            // 1 << 39 (53 - 13)
            if (sticky) dx.x |= (1L << (53L - numPrec));
            mpfr_clear(r);
            return dx.d;
            
        } else {
            // Value is negative
            
            // Take care of denormalized values:
            // 1000 0000 0000 0000 0000 0000 0000 0000 00 = -0
            // 1000 0000 0000 0000 0000 0000 0000 0000 01 = -2^-151          1
            
            // From here on out, we can use a pattern.
            // 1000 0000 0000 0000 0000 0000 0000 0000 10 = -2^-150          2
            // 1000 0000 0000 0000 0000 0000 0000 0000 11 = -1.5 * 2^-150    2
            // 1000 0000 0000 0000 0000 0000 0000 0001 00 = -2^-149          3
            // 1000 0000 0000 0000 0000 0000 0000 0010 00 = -2^-148          4
            // 1000 0000 0000 0000 0000 0000 0000 0100 00 = -2^-147          5
            // 1000 0000 0000 0000 0000 0000 0000 1000 00 = -2^-146          6
            // 1000 0000 0000 0000 0000 0000 0001 0000 00 = -2^-145          7
            // 1000 0000 0000 0000 0000 0000 0010 0000 00 = -2^-144          8
            // 1000 0000 0000 0000 0000 0000 0100 0000 00 = -2^-143          9
            // 1000 0000 0000 0000 0000 0000 1000 0000 00 = -2^-142          10
            // 1000 0000 0000 0000 0000 0001 0000 0000 00 = -2^-141          11
            // 1000 0000 0000 0000 0000 0010 0000 0000 00 = -2^-140          12
            // 1000 0000 0000 0000 0000 0100 0000 0000 00 = -2^-139          13
            // 1000 0000 0000 0000 0000 1000 0000 0000 00 = -2^-138          14
            // 1000 0000 0000 0000 0001 0000 0000 0000 00 = -2^-137          15
            // 1000 0000 0000 0000 0010 0000 0000 0000 00 = -2^-136          16
            // 1000 0000 0000 0000 0100 0000 0000 0000 00 = -2^-135          17
            // 1000 0000 0000 0000 1000 0000 0000 0000 00 = -2^-134          18
            // 1000 0000 0000 0001 0000 0000 0000 0000 00 = -2^-133          19
            // 1000 0000 0000 0010 0000 0000 0000 0000 00 = -2^-132          20
            // 1000 0000 0000 0100 0000 0000 0000 0000 00 = -2^-131          21
            // 1000 0000 0000 1000 0000 0000 0000 0000 00 = -2^-130          22
            // 1000 0000 0001 0000 0000 0000 0000 0000 00 = -2^-129          23
            // 1000 0000 0010 0000 0000 0000 0000 0000 00 = -2^-128          24
            // 1000 0000 0100 0000 0000 0000 0000 0000 00 = -2^-127          25
            // # precision = (152 + exp)
            // First normalized value:
            // 1000 0000 1000 0000 0000 0000 0000 0000 00 = -2^-126          26
            
            
            // 1000 0000 0000 0000 0000 0000 0000 0000 01 = -2^-151
            // 1000 0000 0000 0000 0000 0000 0000 0000 10 = -2^-150
            // if _mval > -2^-150, then round to -2^-151.
            ldexpVal = ldexp(-1.0, -150);
            retVal = ldexp(-1.0, -151);
            cmp = mpfr_cmp_d(mval, ldexpVal);
            if (cmp > 0) return retVal;
            
            // 1000 0000 0000 0000 0000 0000 0000 0000 10 = -2^-150
            // if _mval == -2^-150 and sticky == 0, then it's exact: -2^-150
            retVal = ldexpVal;
            if (cmp == 0 && sticky == 0) return retVal;
            
            // 1000 0000 0000 0000 0000 0000 0000 0000 10 = -2^-150
            // 1000 0000 0000 0000 0000 0000 0000 0000 11 = -1.5 * 2^-150
            // 1000 0000 0000 0000 0000 0000 0000 0001 00 = -2^-149
            // if _mval > -2^-149, then round to -1.5 * 2^-150
            ldexpVal = ldexp(-1.0, -149);
            retVal = ldexp(-1.5, -150);
            cmp = mpfr_cmp_d(mval, ldexpVal);
            if (cmp > 0) return retVal;
            
            // 1000 0000 0000 0000 0000 0000 0000 0001 00 = -2^-149
            // if _mval == -2^-149 and sticky == 0, then it's exact: -2^-149
            retVal = ldexpVal;
            if (cmp == 0 && sticky == 0) return retVal;
            
            // From here on out, we have at least 3 precision bits, so we can work
            // with the normal logic
            
            // 1. Find out how many precision bits we will have
            long exp;
            double fr = mpfr_get_d_2exp(&exp, mval, MPFR_RNDZ);
            long origExp = exp;
            fr *= 2;
            exp--;
            
            // 2. Round _mval to a value with (numPrec - 1) precision bits using
            // round towards zero. Keep track of sticky bit
            long numPrec = 152 + exp;
            mpfr_t r;
            mpfr_init2(r, numPrec - 1);
            status = mpfr_set(r, mval, MPFR_RNDZ);
            if (status != 0) sticky |= 0x1;
            
            // 3. Convert "r" to double. It should be exactly representable.
            doubleX dx;
            dx.d = mpfr_get_d(r, MPFR_RNDZ);
            
            // 4. Add "sticky bit" to the numPrec-th bit
            // Let's say we have 13 precision bits.
	    // Then since it's denormal in 34-bit representation, double rep has 12 mantissa bits.
            // seee eeee eeee xxxx xxxx xxxx ____ ____ ____ ____ ____ ____ ____ ____ ____ ____
            // 1 << 39 (53 - 13)
            if (sticky) dx.x |= (1L << (53L - numPrec));
            mpfr_clear(r);
            return dx.d;
        }
    }
    
    // Assuming not special / denormal
    // 1. Round _mval to a value with 25 precision bits using
    // round towards zero. Keep track of sticky bit
    mpfr_t r;
    mpfr_init2(r, 25);
    status = mpfr_set(r, mval, MPFR_RNDZ);
    if (status != 0) sticky |= 0x1;
    
    // 3. Convert "r" to double. It should be exactly representable.
    doubleX dx;
    dx.d = mpfr_get_d(r, MPFR_RNDZ);
    
    // 4. Add "sticky bit" to the numPrec-th bit
    // Let's say we have 13 precision bits.
    // seee eeee eeee xxxx xxxx xxxx x___ ____ ____ ____ ____ ____ ____ ____ ____ ____
    // 1 << 39 (52 - 13)
    if (sticky) dx.x |= (1L << 27L);
    mpfr_clear(r);
    return dx.d;
}
