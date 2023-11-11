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
  
  return f * rlibm_OneByF[FIndex];
  
}


bool compute_special_case(float x, double& res){

  float_x inp = {.f = x};
  uint32_t ux = inp.x;
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

  /* power of 2 */
  if(__builtin_expect(!m, 0)) {
      res= exp;
      return true;
    }
  return false;
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

  return yp +  rlibm_log2F[FIndex] + exp;
  
}


			       
int main(int argc, char** argv){
  mpfr_init2(mval, 200);

  if(argc != 3){
    printf("Usage: %s <Name of the Interval File> <Oracle File>\n", argv[0]);
    exit(0);
  }

  create_interval_file(argv[1], argv[2], 0x0llu, 0x100000000llu);
  mpfr_clear(mval);
}
