#include<stdio.h>
#include "rlibm.h"
#include "IntervalGen.h"
#include <math.h>
#include "exp2.h"
#include <cassert>
#include <cstdio>



extern mpfr_t mval;

double OutputCompensation(float, double);

void create_interval_file(char*, char*,
			  unsigned long long,
			  unsigned long long);


bool compute_special_case(float x, double& res) {

    float_x fx;
    fx.f = x;

    // 1.  If x == 0x0                      -> 1.0
    // 2.  If 0 < x <= 864026623            -> 1.0 + 1 ulp
    // 3.  If 864026623 < x <= 872415231    -> 1.0 + 3 ulp
    // 4.  If 1118925336 <= x < 0x7F800000  -> MAXVAL
    // 5.  If x == 0x7F800000               -> Infinity
    // 6.  If 0x7F800000 < x < 0x80000000   -> NaN

    // 7.  If x == 0x80000000               -> 1.0
    // 8.  If 0x80000000 < x <= 3003121664  -> 1.0 - 1 ulp
    // 9.  If 3003121664 < x <= 3011510272  -> 1.0 - 3 ulp
    // 10. If 3268407733 <= x < 0xFF800000  -> 2^-151
    // 11. If x == 0xFF800000               -> 0.0
    // 12. If x > 0xFF800000                -> NaN


    if ((fx.x & 0x7FFFFFFF) == 0) {
      res = 1.0;
      return true;
    }

    if (fx.x <= 872415231) {
      if (fx.x <= 864026623) {
	res = 1.0000000298023223876953125;
	return true;
      }
      res = 1.0000000894069671630859375;
      return true;
    }

    if (1118925336 <= fx.x && fx.x <= 3011510272) {
      if (fx.x < 0x80000000) {
	// positive counterpart
	if (fx.x < 0x7F800000) {
	  res = 3.40282361850336062550457001444955389952e+38;
	  return true;
	}
	if (fx.x == 0x7F800000) {
	  res = 1.0 / 0.0;
	  return true;
	}
	res = 0.0/0.0;
	return true;
      }
      
      // negative counterpart
      if (fx.x <= 3003121664) {
	res = 0.99999998509883880615234375;
	return true;
      }

      res = 0.99999995529651641845703125;
      return true;
    }

    if (fx.x >= 3268407733) {
      if (fx.x == 0xFF800000) {
        res = 0.0;
	return true;
      }
      if (fx.x < 0xFF800000) {
	res = ldexp(1.0, -151);
	return true;
      }
      res = 0.0/0.0;
      return true;
    }

    return false;

}

double RangeReduction(float x) {
    double xp = x * 92.332482616893656768297660164535045623779296875;
    int N = (int)xp;
    
    return x - N *
    0.01083042469624914509729318723429969395510852336883544921875;

}
    
double OutputCompensation(float x, double yp) {
    double xp = x * 92.332482616893656768297660164535045623779296875;
    int N = (int)xp;
    int N2 = N % 64;
    if (N2 < 0) N2 += 64;
    int N1 = N - N2;
    int M = N1 / 64;
    
    return yp * ldexp(exp2JBy64[N2], M);
}

void GuessInitialLbUb(double x,
		      double roundingLb, double roundingUb,
		      double xp, double& lb, double& ub) {
  // Take a guess of yp that will end up in roundingLb, roundingUb.
  double_x tempYp;
  tempYp.d = exp(xp);
  double tempY = OutputCompensation(x, tempYp.d);
    
  if (tempY < roundingLb) {
    // if tempY < roundingLb, then keep increasing tempYp until tempY is
    // >= roundingLb.
    do {
      if (tempYp.d >= 0.0) tempYp.x++;
      else tempYp.x--;
        tempY = OutputCompensation(x, tempYp.d);
    } while (tempY < roundingLb);
        
    // Then, it had better be that roundingLb <= tempY <= roundingUb.
    if (tempY > roundingUb) {
      printf("Error during GuessInitialLbUb: lb > ub.\n");
      printf("x = %.100e\n", x);
      exit(0);
    }
    lb = tempYp.d;
    ub = tempYp.d;
    return;
  }
    
  if (tempY > roundingUb) {
    // if tempY > roundingUb, then keep decreasing tempYp until tempY is
    // <= roundingUb.
    do {
      if (tempYp.d >= 0.0) tempYp.x--;
      else tempYp.x++;
      tempY = OutputCompensation(x, tempYp.d);
    } while (tempY > roundingUb);
        
    // Then, it had better be that roundingLb <= tempY <= roundingUb.
    if (tempY < roundingLb) {
      printf("Error during GuessInitialLbUb: lb > ub.\n");
      printf("x = %.100e\n", x);
      exit(0);
    }
    lb = tempYp.d;
    ub = tempYp.d;
    return;
  }
    
  lb = tempYp.d;
  ub = tempYp.d;
  return;
}

void SpecCaseRedInt(float x,
		    double glb, bool& blb, double& slb,
		    double gub, bool& bub, double& sub) {
    blb = false;
    bub = false;
    return;
}

int main(int argc, char** argv) {
    mpfr_init2(mval, 200);
    
    if (argc != 3) {
        printf("Usage: %s <Name of the Interval File> <Oracle File>\n", argv[0]);
        exit(0);
    }

    create_interval_file(argv[1], argv[2], 0x0llu, 0x100000000llu);
    mpfr_clear(mval);

    return 0;
}
