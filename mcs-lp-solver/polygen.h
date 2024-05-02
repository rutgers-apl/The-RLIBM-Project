#pragma once

#define SOPLEX_WITH_GMP
#include <math.h>
#include <set>
#include<stdio.h>
#include<stdlib.h>
#include <soplex.h>
#include<fstream>
#include<memory>
#include<random>

#define MAX_TRIES 100
#define VIOLATE_THRESHOLD 0
#define SAMPLE_MATCH_THRESHOLD 20
#define MAX_ITERATIONS 1200

#include "maxfs_common.h"

using namespace soplex;
using namespace std;

typedef struct {
  double x;         /* original input */
  double lb;        /* lower bound */
  double ub;        /* upper bound */
  double w;         /* weight */
  double u;         /* uniform random value */
} interval_data;

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
  double key;
  size_t index;
} sample_info;


typedef union {
  double d;
  uint64_t x;
} double_x;

double rlibm_poly_evaluation(double, polynomial*);
