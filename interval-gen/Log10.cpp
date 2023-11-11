#include<stdio.h>
#include "rlibm.h"
#include "IntervalGen.h"
#include <math.h>
#include "rlibm_log.h"
#include <cassert>
#include <cstdio>


extern mpfr_t mval;

#define LOG102HIGH 0x1.34413509f79fep-2
#define LOG102LOW  0x1.e623e2566b02ep-55



double OutputCompensation(float, double);

void create_interval_file(char*, char*,
			  unsigned long long,
			  unsigned long long);


bool compute_special_case(float x, double& res) {

  float_x inp = {.f = x};
  uint32_t ux = inp.x;
  
  switch (ux) {
  case 0x3f800000 :
    res = 0.0;
    return true;
  case 0x41200000 :
    res = 1.0;
    return true;
  case 0x42c80000 :
    res = 2.0;
    return true;
  case 0x447a0000 :
    res = 3.0;
    return true;
  case 0x461c4000 :
    res = 4.0;
    return true;
  case 0x47c35000 :
    res = 5.0;
    return true;
  case 0x49742400:
    res = 6.0;
    return true;
  case 0x4b189680 :
    res = 7.0;
    return true;
  case 0x4cbebc20 :
    res = 8.0;
    return true;
  case 0x4e6e6b28 :
    res = 9.0;
    return true;
  case 0x501502f9 :
    res = 10.0;
    return true;
  }

  
  uint64_t m = ux & 0x7FFFFF;
  m = m << 29;
  int exp = (ux >> 23) - 127;
  
  if(__builtin_expect(ux < 0x800000 || ux >= 0x7F800000, 0)){

    if (ux==0||ux==(1u<<31)){
      res = -__builtin_inff(); // +0.0 || -0.0
      return true;
    }

    uint32_t inf_or_nan = ((ux>>23)&0xff) == 0xff, nan = inf_or_nan && (ux<<9);

    if (ux>>31 && !nan) {
     res =  __builtin_nanf("-");
     return true;
    }

    if (inf_or_nan) {
      res = x;
      return true;
    }

    // subnormal
    int nz = __builtin_clzll(m);
    m <<= nz-11;
    m &= ~0ul>>12;
    exp = exp - (nz - 12);
  }

  return false;

  
}

double RangeReduction(float x) {
  
  float_x inp = {.f = x};
  uint32_t ux = inp.x;
  uint64_t m = ux & 0x7FFFFF;
  m = m << 29;
  int exp = (ux >> 23) - 127;
  
  if(__builtin_expect(ux < 0x800000 || ux >= 0x7F800000, 0)){
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

  double_x new_rr= {.d = f * rlibm_OneByF[FIndex]};

  return new_rr.d;

  
}

    
double OutputCompensation(float x, double yp) {

  float_x inp = {.f = x};
  uint32_t ux = inp.x;
  uint64_t m = ux & 0x7FFFFF;
  m = m << 29;
  int exp = (ux >> 23) - 127;
  
  if(__builtin_expect(ux < 0x800000 || ux >= 0x7F800000, 0)){
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
  
  f *= rlibm_OneByF[FIndex];

  /* we are approximating log1p(f) with the polynomial */

  double_x new_oc = {.d =  yp + exp * LOG102LOW + rlibm_log10F[FIndex] + exp * LOG102HIGH};

  return new_oc.d;

}

void GuessInitialLbUb(double x,
		      double roundingLb, double roundingUb,
		      double xp, double& lb, double& ub) {
  // Take a guess of yp that will end up in roundingLb, roundingUb.
  double_x tempYp;
  tempYp.d = log1p(xp) /log(10);
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
