#include <stdio.h>
#include "rlibm.h"
#include "IntervalGen.h"
#include <math.h>
#include "rlibm_log.h"
#include <cassert>
#include <cstdio>

extern mpfr_t mval;

double OutputCompensation(float, double);

void create_interval_file(char*, char*,
			  unsigned long long,
			  unsigned long long);

void GuessInitialLbUb(double x,
		      double roundingLb, double roundingUb,
		      double xp, double& lb, double& ub) {
  // Take a guess of yp that will end up in roundingLb, roundingUb.                                
  double_x tempYp;
  tempYp.d = log1p(xp) / log(2);
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
      printf("x = %a\n", x);
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


double RangeReduction(float x) {
  float_x fix, fit;
  
  int m = 0;
  fix.f = x;
  if (fix.x < 0x800000) {
    fix.f *= pow(2, 23);
    m -= 23;
  }
  m += fix.x >> 23;
  m -= 127;
  fix.x &= 0x007FFFFF;
  fix.x |= 0x3F800000;
  
  fit.x = fix.x & 0x007F0000;
  int FIndex = fit.x >> 16;
  fit.x |= 0x3F800000;
  double F = fit.f;
  
  double f = fix.f - F;
  return f * rlibm_OneByF[FIndex];
}


bool compute_special_case(float x, double& res){

  float_x fx;
  fx.f = x;
  if (x == 0.0) {
    res = -1.0/0.0;
    return true;
  } else if (fx.x == 0x7F800000) {
    res = x;
    return true;
  } else if (fx.x > 0x7F800000) {
    fx.x = 0x7FFFFFFF;
    res = fx.f;
    return true;
  }  
  int exp;
  float remainder = frexpf(fx.f, &exp);
  
  if (remainder == 0.5f || remainder == -0.5f) {
    res = (exp - 1);
    return true;
  }
  
  return false;
}

double OutputCompensation(float x, double yp) {
  float_x fix, fit;
  
  int m = 0;
  fix.f = x;
  if (fix.x < 0x800000) {
    fix.f *= pow(2, 23);
    m -= 23;
  }
  m += fix.x >> 23;
  m -= 127;
  fix.x &= 0x007FFFFF;
  fix.x |= 0x3F800000;
  
  fit.x = fix.x & 0x007F0000;
  int FIndex = fit.x >> 16;
  
  return yp + rlibm_log2F[FIndex] + m;
}


			       
int main(int argc, char** argv){
  mpfr_init2(mval, 200);

  if(argc != 3){
    printf("Usage: %s <Name of Interval File> <Oracle File>\n", argv[0]);
    exit(0);
  }

  create_interval_file(argv[1], argv[2], 0x0llu, 0x100000000llu);
  mpfr_clear(mval);
}
