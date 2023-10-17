#pragma once

#define SOPLEX_WITH_GMP
#include <math.h>
#include<stdio.h>
#include<stdlib.h>
#include <soplex.h>
#include<fstream>
#include<memory>
#include<random>

#define MAX_TRIES 100
#define VIOLATE_THRESHOLD 10
#define SAMPLE_MATCH_THRESHOLD 20
#define MAX_ITERATIONS 600





#ifdef EXIT_ON_THRESHOLD
const int RLIBM_EXIT_ON_THRESHOLD = 1;
#else
const int RLIBM_EXIT_ON_THRESHOLD = 0;
#endif

using namespace soplex;
using namespace std;


typedef struct {
  double x;         /* original input */
  double lb;        /* lower bound */
  double ub;        /* upper bound */
  double w;         /* weight */
  double u;         /* uniform random value */
} interval_data;


/* progressive interval data */
typedef struct {
  double x;         /* original input */
  double lb;        /* lower bound */
  double ub;        /* upper bound */
  double w;         /* weight */
  double u;         /* uniform random value */
  int rep;          /* representation bfloat16 (0) , tensorfloat32 (1), float32 (2) */
} pinterval_data;


typedef struct{
  double x;
  double lb;
  double ub;
  double orig_lb;  /* remember lb before narrowing */
  double orig_ub;  /* remember ub before narrowing */
  double w;
  double u;
  double k;        /* key computed as 1/u^w */

} sample_data;


typedef struct{
  double x;
  double lb;
  double ub;
  double orig_lb;  /* remember lb before narrowing */
  double orig_ub;  /* remember ub before narrowing */
  double w;
  double u;
  double k;        /* key computed as 1/u^w */
  int rep;          /* representation bfloat16 (0) , tensorfloat32 (1), float32 (2) */
} psample_data;


typedef struct{
  double key;
  size_t index;
} sample_info;

typedef struct {
  int termsize;
  int* power;
  double* coeffs;
} polynomial;

typedef union {
  double d;
  uint64_t x;
} double_x;

