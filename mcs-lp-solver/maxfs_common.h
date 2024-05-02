#pragma once


#include <math.h>
#include <cassert>
#include <cstring>

typedef enum {
  NONE=0,
  LOG2=1,
  LOG=2,
  LOG10=3,
  EXP2=4,
  EXP=5,
  EXP10=6,
  SINH_SINH=7,
  SINH_COSH=8,
  COSH_SINH=9,
  COSH_COSH=10,
  SINPI_SINPI=11,
  SINPI_COSPI=12,
  COSPI_SINPI=13,
  COSPI_COSPI=14,
  SIN_SIN=15,
  SIN_COS=16,
  COS_SIN=17,
  COS_COS=18,  
} my_function;


typedef struct {
  int termsize;
  int* power;
  double* coeffs;
} polynomial;

double rlibm_poly_evaluation(double, polynomial*);

extern my_function function_to_process;

void set_function_process(char**);

