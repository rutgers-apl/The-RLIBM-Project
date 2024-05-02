#define SOPLEX_WITH_GMP
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <soplex.h>
#include <chrono>
#include <fstream>
#include <memory>
#include <random>
#include <unordered_set>

#include "maxfs_common.h"

#define MAX_TRIES 200
int VIOLATE_THRESHOLD=11;
#define SAMPLE_MATCH_THRESHOLD 20
#define MAX_ITERATIONS 1200

#define EXIT_ON_THRESHOLD 1

#ifdef EXIT_ON_THRESHOLD
const int RLIBM_EXIT_ON_THRESHOLD = 1;
#else
const int RLIBM_EXIT_ON_THRESHOLD = 0;
#endif

my_function function_to_process = NONE;

using namespace soplex;
using namespace std;
using namespace chrono;

size_t prev_successful_degree = 0;

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


polynomial* rlibm_solve_with_soplex(sample_data* sintervals, size_t ssize,
                                    sample_data* sp_intervals, size_t sp_size,
                                    int* power, int termsize){
  
  double* u = (double*) calloc(sp_size, sizeof(double));      // indicator variable for upper bound
  double* w = (double*) calloc(sp_size, sizeof(double));      // indicator variable for lower bound
  double* s = (double*) calloc(sp_size, sizeof(double));      // slack variable for upper bound
  double* r = (double*) calloc(sp_size, sizeof(double));      // slack variable for lower bound

  // initialize indicator/slack variables
  for(unsigned k = 0; k < sp_size; ++k){
      u[k] = 1.;
      w[k] = 1.;
      s[k] = 0.;
      r[k] = 0.;
  }
  
  SoPlex mysoplex;
  mysoplex.setBoolParam(SoPlex::RATFACJUMP,true);
  mysoplex.setIntParam(SoPlex::SOLVEMODE,2);
  mysoplex.setIntParam(SoPlex::CHECKMODE,2);
  mysoplex.setIntParam(SoPlex::SYNCMODE,1);
  mysoplex.setIntParam(SoPlex::READMODE,1);
  mysoplex.setRealParam(SoPlex::FEASTOL,0.0);
  mysoplex.setRealParam(SoPlex::OPTTOL,0.0);
  mysoplex.setRealParam(SoPlex::EPSILON_ZERO,0.0);
  mysoplex.setRealParam(SoPlex::EPSILON_FACTORIZATION,0.0);
  mysoplex.setRealParam(SoPlex::EPSILON_UPDATE,0.0);
  mysoplex.setRealParam(SoPlex::EPSILON_PIVOT,0.0);
  mysoplex.setIntParam(SoPlex::VERBOSITY,0);
  mysoplex.setRealParam(SoPlex::TIMELIMIT,5*60);

  // set objective sense
  mysoplex.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MINIMIZE);
  
  /* we first add objective variables */
  DSVectorRational dummycol(0);

  // objective function for slack variables
  for(unsigned long k = 0; k < sp_size; ++k){
      auto column1 = LPColRational(1.0, dummycol, infinity, 0.);
      mysoplex.addColRational(column1);
      auto column2 = LPColRational(1.0, dummycol, infinity, 0.);
      mysoplex.addColRational(column2);
  }

  // objective function for polynomial coefficients
  for(int i = 0; i < termsize; i++){
      Rational coeff(0.0);
      for(unsigned long k = 0; k < sp_size; ++k){
          Rational xValR(sp_intervals[k].x);

          Rational toAdd(1.0);
          for(int t = 0; t < power[i]; ++t) toAdd*=xValR;

          coeff += Rational(w[k] - u[k])*toAdd;
      }
      auto column = LPColRational(coeff, dummycol, infinity, -infinity);
      mysoplex.addColRational(column);
  }

  // add linear constraints for special cases
  for(unsigned long k = 0; k < sp_size; ++k){
      DSVectorRational row1(2*sp_size + termsize);
      DSVectorRational row2(2*sp_size + termsize);
      Rational xValR(sp_intervals[k].x);

      row1.add(2*k, Rational(1.0));       // s
      row2.add(2*k+1, Rational(1.0));     // r

      for(int j=0; j < termsize;j++){
          Rational toAdd(1.0);
          for(int t = 0; t < power[j]; t++) toAdd*=xValR;

          row1.add(2*sp_size + j, -toAdd);
          row2.add(2*sp_size + j, toAdd);
      }

      double lbnd = sp_intervals[k].lb;
      double ubnd = sp_intervals[k].ub;
      mysoplex.addRowRational(LPRowRational(-ubnd,row1,infinity));
      mysoplex.addRowRational(LPRowRational(lbnd,row2,infinity));
  }

  // add linear constraints for basis intervals
  for(int i = 0; i < ssize; i++){
    Rational xValR(sintervals[i].x);
    DSVectorRational row1(2*sp_size + termsize);
    
    for(int j=0; j<termsize;j++){
      Rational toAdd(1.0);
      for(int k=0;k<power[j];k++) toAdd*=xValR;

      row1.add(2*sp_size + j,toAdd);
    }
        
    // LPRow: low bound, row, upper bound
    double lbnd= sintervals[i].lb;
    double ubnd= sintervals[i].ub;
    mysoplex.addRowRational(LPRowRational(lbnd,row1,ubnd));
  }

  /* solve LP */
  SPxSolver::Status stat;
  stat=mysoplex.optimize();
  
  /* get solution */
  if(stat==SPxSolver::OPTIMAL){
    DVectorRational prim(2*sp_size + termsize);
    mysoplex.getPrimalRational(prim);

    /* generate the polynomial as a plain structure */
    polynomial* p = (polynomial *) calloc(1, sizeof(polynomial));
    p->termsize = termsize;
    p->power = power;
    p->coeffs = (double *) calloc(termsize, sizeof(double));

    for(int i = 0; i < termsize; i++)
      p->coeffs[i] = mpq_get_d(*(prim[2*sp_size + i].getMpqPtr_w()));
    
    return p;
  }
  else if(stat==SPxSolver::UNBOUNDED){

    polynomial* p = (polynomial *) calloc(1, sizeof(polynomial));
    p->termsize = termsize;
    p->power = power;
    p->coeffs = (double *) calloc(termsize, sizeof(double));
    
    for(int i=0;i<termsize;i++)
      p->coeffs[i] = 0.0;
    
    return p;
  }

  free(u);
  free(w);
  free(s);
  free(r);
  
  return nullptr;
}

void check_sorted(sample_info* sampled_indices, size_t ssize){
  double min= sampled_indices[0].key;

  for(size_t i = 0; i< ssize; i++){
    assert ( min <= sampled_indices[i].key);
    min = sampled_indices[i].key;
  }
  
}

# if 0

//moved to polyeval
double rlibm_poly_evaluation(double x, polynomial* poly){

#if 1   
    double  C1 = poly->coeffs[0];
    double  C2 = poly->coeffs[1];
    double  C3 = poly->coeffs[2];
    double  C4 = poly->coeffs[3];
    double  C5 = poly->coeffs[4];
    double  C6 = poly->coeffs[5];

      double xsquare = x * x;
  double xcube = x * xsquare;
  double temp1 = fma(C2, x, C1);
  //  double temp1 = C2 * x +  C1;
  double temp2 = fma(C3, xsquare, temp1);
  
  //  double temp2 = C3 * xsquare + temp1;

  double temp3 = fma(C5, x, C4);
  
  //  double temp3 = C5 * x + C4;

  double temp4 = fma(C6, xsquare, temp3);
  //  double temp4 = C6 * xsquare + temp3;

  double temp5 = fma(temp4, xcube, temp2);
  //  double temp5 = temp4 *xcube + temp2;

  return x * temp5;


#endif

#if 0
  // for even functions
  if(poly->termsize == 1){
      return poly->coeffs[0];
  } else if(poly->termsize == 2){
    return fma(x*x, poly->coeffs[1], poly->coeffs[0]);
  } else if(poly->termsize == 3){
    double xsquare = x*x;
    double temp = fma(xsquare, poly->coeffs[2], poly->coeffs[1]);
    return fma(xsquare, temp, poly->coeffs[0]);
  } else if(poly->termsize == 4){
    double xsquare = x*x;
    double temp1 = fma(xsquare, poly->coeffs[1], poly->coeffs[0]);
    double temp2 = fma(xsquare, poly->coeffs[3], poly->coeffs[2]);
    return fma(xsquare*xsquare, temp2, temp1);
  } else if(poly->termsize == 5){
    double x2 = x*x, x4 = x2*x2;
    double temp1 = fma(x2, poly->coeffs[1], poly->coeffs[0]);
    double temp2 = fma(x2, poly->coeffs[3], poly->coeffs[2]);
    double temp3 = fma(x4, poly->coeffs[4], temp2);
    return fma(x4, temp3, temp1);
  } else if(poly->termsize == 6){
    double x2 = x*x, x4 = x2*x2;
    double temp1 = fma(x2, poly->coeffs[1], poly->coeffs[0]);
    double temp2 = fma(x2, poly->coeffs[3], poly->coeffs[2]);
    double temp3 = fma(x2, poly->coeffs[5], poly->coeffs[4]);
    double temp4 = fma(x4, temp3, temp2);
    return fma(x4, temp4, temp1);
  } else if(poly->termsize == 7){
    double x2 = x*x, x4 = x2*x2;
    double temp1 = fma(x2, poly->coeffs[1], poly->coeffs[0]);
    double temp2 = fma(x2, poly->coeffs[3], poly->coeffs[2]);
    double temp3 = fma(x2, poly->coeffs[5], poly->coeffs[4]);
    double temp4 = fma(x4, poly->coeffs[6], temp3);
    double temp5 = fma(x4, temp4, temp2);
    return fma(x4, temp5, temp1);
  } else if(poly->termsize == 8){
    double x2 = x*x, x4 = x2*x2;
    double temp1 = fma(x2, poly->coeffs[1], poly->coeffs[0]);
    double temp2 = fma(x2, poly->coeffs[3], poly->coeffs[2]);
    double temp3 = fma(x2, poly->coeffs[5], poly->coeffs[4]);
    double temp4 = fma(x2, poly->coeffs[7], poly->coeffs[6]);
    double temp5 = fma(x4, temp4, temp3);
    double temp6 = fma(x4, temp5, temp2);
    return fma(x4, temp6, temp1);
  }
#endif

#if 0
  // for odd functions
  if(poly->termsize == 1){
      return x*poly->coeffs[0];
  } else if(poly->termsize == 2){
    return x*fma(x*x, poly->coeffs[1], poly->coeffs[0]);
  } else if(poly->termsize == 3){
    double xsquare = x*x;
    double temp = fma(xsquare, poly->coeffs[2], poly->coeffs[1]);
    return x*fma(xsquare, temp, poly->coeffs[0]);
  } else if(poly->termsize == 4){
    double xsquare = x*x;
    double temp1 = fma(xsquare, poly->coeffs[1], poly->coeffs[0]);
    double temp2 = fma(xsquare, poly->coeffs[3], poly->coeffs[2]);
    return x*fma(xsquare*xsquare, temp2, temp1);
  } else if(poly->termsize == 5){
    double x2 = x*x, x4 = x2*x2;
    double temp1 = fma(x2, poly->coeffs[1], poly->coeffs[0]);
    double temp2 = fma(x2, poly->coeffs[3], poly->coeffs[2]);
    double temp3 = fma(x4, poly->coeffs[4], temp2);
    return x*fma(x4, temp3, temp1);
  } else if(poly->termsize == 6){
    double x2 = x*x, x4 = x2*x2;
    double temp1 = fma(x2, poly->coeffs[1], poly->coeffs[0]);
    double temp2 = fma(x2, poly->coeffs[3], poly->coeffs[2]);
    double temp3 = fma(x2, poly->coeffs[5], poly->coeffs[4]);
    double temp4 = fma(x4, temp3, temp2);
    return x*fma(x4, temp4, temp1);
  } else if(poly->termsize == 7){
    double x2 = x*x, x4 = x2*x2;
    double temp1 = fma(x2, poly->coeffs[1], poly->coeffs[0]);
    double temp2 = fma(x2, poly->coeffs[3], poly->coeffs[2]);
    double temp3 = fma(x2, poly->coeffs[5], poly->coeffs[4]);
    double temp4 = fma(x4, poly->coeffs[6], temp3);
    double temp5 = fma(x4, temp4, temp2);
    return x*fma(x4, temp5, temp1);
  } else if(poly->termsize == 8){
    double x2 = x*x, x4 = x2*x2;
    double temp1 = fma(x2, poly->coeffs[1], poly->coeffs[0]);
    double temp2 = fma(x2, poly->coeffs[3], poly->coeffs[2]);
    double temp3 = fma(x2, poly->coeffs[5], poly->coeffs[4]);
    double temp4 = fma(x2, poly->coeffs[7], poly->coeffs[6]);
    double temp5 = fma(x4, temp4, temp3);
    double temp6 = fma(x4, temp5, temp2);
    return x*fma(x4, temp6, temp1);
  }
#endif

#if 0
  // for log
  if(poly->termsize == 1){
      return x*poly->coeffs[0];
  } else if(poly->termsize == 2){
    return x*fma(x, poly->coeffs[1], poly->coeffs[0]);
  } else if(poly->termsize == 3){
    double temp = fma(x, poly->coeffs[2], poly->coeffs[1]);
    return x*fma(x, temp, poly->coeffs[0]);
  } else if(poly->termsize == 4){
    double x2 = x*x;
    double temp1 = fma(x, poly->coeffs[1], poly->coeffs[0]);
    double temp2 = fma(x, poly->coeffs[3], poly->coeffs[2]);
    return x*fma(x2, temp2, temp1);
  } else if(poly->termsize == 5){
    double x2 = x*x;
    double temp1 = fma(x, poly->coeffs[1], poly->coeffs[0]);
    double temp2 = fma(x, poly->coeffs[3], poly->coeffs[2]);
    double temp3 = fma(x2, poly->coeffs[4], temp2);
    return x*fma(x2, temp3, temp1);
  } else if(poly->termsize == 6){
    double x2 = x*x;
    double temp1 = fma(x, poly->coeffs[1], poly->coeffs[0]);
    double temp2 = fma(x, poly->coeffs[3], poly->coeffs[2]);
    double temp3 = fma(x, poly->coeffs[5], poly->coeffs[4]);
    double temp4 = fma(x2, temp3, temp2);
    return x*fma(x2, temp4, temp1);
  } else if(poly->termsize == 7){
    double x2 = x*x;
    double temp1 = fma(x, poly->coeffs[1], poly->coeffs[0]);
    double temp2 = fma(x, poly->coeffs[3], poly->coeffs[2]);
    double temp3 = fma(x, poly->coeffs[5], poly->coeffs[4]);
    double temp4 = fma(x2, poly->coeffs[6], temp3);
    double temp5 = fma(x2, temp4, temp2);
    return x*fma(x2, temp5, temp1);
  } else if(poly->termsize == 8){
    double x2 = x*x;
    double temp1 = fma(x, poly->coeffs[1], poly->coeffs[0]);
    double temp2 = fma(x, poly->coeffs[3], poly->coeffs[2]);
    double temp3 = fma(x, poly->coeffs[5], poly->coeffs[4]);
    double temp4 = fma(x, poly->coeffs[6], poly->coeffs[7]);
    double temp5 = fma(x2, temp4, temp3);
    double temp6 = fma(x2, temp5, temp2);
    return x*fma(x2, temp6, temp1);
  } else if(poly->termsize == 9){
    double x2 = x*x;
    double temp1 = fma(x, poly->coeffs[1], poly->coeffs[0]);
    double temp2 = fma(x, poly->coeffs[3], poly->coeffs[2]);
    double temp3 = fma(x, poly->coeffs[5], poly->coeffs[4]);
    double temp4 = fma(x, poly->coeffs[6], poly->coeffs[7]);
    double temp5 = fma(x2, poly->coeffs[8], temp4);
    double temp6 = fma(x2, temp5, temp3);
    double temp7 = fma(x2, temp6, temp2);
    return x*fma(x2, temp7, temp1);
  } else if(poly->termsize == 10){
    double x2 = x*x;
    double temp1 = fma(x, poly->coeffs[1], poly->coeffs[0]);
    double temp2 = fma(x, poly->coeffs[3], poly->coeffs[2]);
    double temp3 = fma(x, poly->coeffs[5], poly->coeffs[4]);
    double temp4 = fma(x, poly->coeffs[6], poly->coeffs[7]);
    double temp5 = fma(x, poly->coeffs[8], poly->coeffs[9]);
    double temp6 = fma(x2, temp5, temp4);
    double temp7 = fma(x2, temp6, temp3);
    double temp8 = fma(x2, temp7, temp2);
    return x*fma(x2, temp8, temp1);
  }
#endif

#if 0
  // for exp
  if(poly->termsize == 1){
      return poly->coeffs[0];
  } else if(poly->termsize == 2){
    return fma(x, poly->coeffs[1], poly->coeffs[0]);
  } else if(poly->termsize == 3){
    double temp = fma(x, poly->coeffs[2], poly->coeffs[1]);
    return fma(x, temp, poly->coeffs[0]);
  } else if(poly->termsize == 4){
    double x2 = x*x;
    double temp1 = fma(x, poly->coeffs[1], poly->coeffs[0]);
    double temp2 = fma(x, poly->coeffs[3], poly->coeffs[2]);
    return fma(x2, temp2, temp1);
  } else if(poly->termsize == 5){
    double x2 = x*x;
    double temp1 = fma(x, poly->coeffs[1], poly->coeffs[0]);
    double temp2 = fma(x, poly->coeffs[3], poly->coeffs[2]);
    double temp3 = fma(x2, poly->coeffs[4], temp2);
    return fma(x2, temp3, temp1);
  } else if(poly->termsize == 6){
    double x2 = x*x;
    double temp1 = fma(x, poly->coeffs[1], poly->coeffs[0]);
    double temp2 = fma(x, poly->coeffs[3], poly->coeffs[2]);
    double temp3 = fma(x, poly->coeffs[5], poly->coeffs[4]);
    double temp4 = fma(x2, temp3, temp2);
    return fma(x2, temp4, temp1);
  } else if(poly->termsize == 7){
    double x2 = x*x;
    double temp1 = fma(x, poly->coeffs[1], poly->coeffs[0]);
    double temp2 = fma(x, poly->coeffs[3], poly->coeffs[2]);
    double temp3 = fma(x, poly->coeffs[5], poly->coeffs[4]);
    double temp4 = fma(x2, poly->coeffs[6], temp3);
    double temp5 = fma(x2, temp4, temp2);
    return fma(x2, temp5, temp1);
  } else if(poly->termsize == 8){
    double x2 = x*x;
    double temp1 = fma(x, poly->coeffs[1], poly->coeffs[0]);
    double temp2 = fma(x, poly->coeffs[3], poly->coeffs[2]);
    double temp3 = fma(x, poly->coeffs[5], poly->coeffs[4]);
    double temp4 = fma(x, poly->coeffs[6], poly->coeffs[7]);
    double temp5 = fma(x2, temp4, temp3);
    double temp6 = fma(x2, temp5, temp2);
    return fma(x2, temp6, temp1);
  } else if(poly->termsize == 9){
    double x2 = x*x;
    double temp1 = fma(x, poly->coeffs[1], poly->coeffs[0]);
    double temp2 = fma(x, poly->coeffs[3], poly->coeffs[2]);
    double temp3 = fma(x, poly->coeffs[5], poly->coeffs[4]);
    double temp4 = fma(x, poly->coeffs[6], poly->coeffs[7]);
    double temp5 = fma(x2, poly->coeffs[8], temp4);
    double temp6 = fma(x2, temp5, temp3);
    double temp7 = fma(x2, temp6, temp2);
    return fma(x2, temp7, temp1);
  } else if(poly->termsize == 10){
    double x2 = x*x;
    double temp1 = fma(x, poly->coeffs[1], poly->coeffs[0]);
    double temp2 = fma(x, poly->coeffs[3], poly->coeffs[2]);
    double temp3 = fma(x, poly->coeffs[5], poly->coeffs[4]);
    double temp4 = fma(x, poly->coeffs[6], poly->coeffs[7]);
    double temp5 = fma(x, poly->coeffs[8], poly->coeffs[9]);
    double temp6 = fma(x2, temp5, temp4);
    double temp7 = fma(x2, temp6, temp3);
    double temp8 = fma(x2, temp7, temp2);
    return fma(x2, temp8, temp1);
  }
#endif

  double ret_val = 0.0;

  // simulated Horner's method
  for(int i = poly->termsize-1; i> 0; i--){
    ret_val = ret_val + poly->coeffs[i];
    double xmul = 1.0;
    for(int j = 0; j < (poly->power[i] - poly->power[i-1]); j++){
      xmul = xmul * x;
    }
    ret_val = ret_val * xmul;	  
  }
  ret_val = ret_val + poly->coeffs[0];
  
  for(int j = 0; j < poly->power[0]; j++){
    ret_val = ret_val * x;
  }  
  return ret_val;
}
#endif

bool rlibm_validate_and_fix_intervals(sample_data* sintervals,
				      size_t ssize, polynomial* poly){

  bool return_val = true;
  for(size_t i = 0; i < ssize; i++){
    double y = rlibm_poly_evaluation(sintervals[i].x, poly);

    if(y < sintervals[i].orig_lb){
      return_val = false;
      double_x lbx;
      lbx.d = sintervals[i].lb;
      if(lbx.d >= 0) {
	lbx.x = lbx.x + 1;
      }
      else{
	lbx.x = lbx.x - 1 ;
      }
      sintervals[i].lb = lbx.d;
    }
    else if(y > sintervals[i].orig_ub){
      return_val = false;
      double_x ubx;
      ubx.d = sintervals[i].ub;
      if(ubx.d >= 0){
	ubx.x = ubx.x - 1;
      }
      else {
	ubx.x = ubx.x + 1;
      }
      sintervals[i].ub = ubx.d;
    }    
  }
  return return_val;
}

size_t rlibm_validate_and_fix_as_many_intervals(sample_data* sintervals,
				      size_t ssize, polynomial* poly){

  size_t num_intervals = 0;
  for(size_t i = 0; i < ssize; i++){
    double y = rlibm_poly_evaluation(sintervals[i].x, poly);

    if(y < sintervals[i].orig_lb){
      ++num_intervals;
      double_x lbx;
      lbx.d = sintervals[i].lb;
      if(lbx.d >= 0) {
	lbx.x = lbx.x + 1;
      }
      else{
	lbx.x = lbx.x - 1 ;
      }
      sintervals[i].lb = lbx.d;
    }
    else if(y > sintervals[i].orig_ub){
      ++num_intervals;
      double_x ubx;
      ubx.d = sintervals[i].ub;
      if(ubx.d >= 0){
	ubx.x = ubx.x - 1;
      }
      else {
	ubx.x = ubx.x + 1;
      }
      sintervals[i].ub = ubx.d;
    }    
  }
  return num_intervals;
}

// memory leak on the polynomial

polynomial*
rlibm_generate_polynomial(sample_data* sintervals, size_t ssize,
                          sample_data* sp_intervals, size_t sp_size,
			  int* power, int power_size, int max_tries){

  for(int i = power_size-1; i < power_size; i++){
    printf("Trying to generate a polynomial  with %d terms \n", i+1);

    int count = 0;
    while(count < max_tries){
      polynomial* p = rlibm_solve_with_soplex(sintervals, ssize, sp_intervals, sp_size, power, i+1);
      //if(p && rlibm_validate_and_fix_intervals(sintervals, ssize, p) && (rlibm_validate_and_fix_as_many_intervals(sp_intervals, sp_size, p) < sp_size)){
      if(p && rlibm_validate_and_fix_intervals(sintervals, ssize, p) && (rlibm_validate_and_fix_as_many_intervals(sp_intervals, sp_size, p) <= VIOLATE_THRESHOLD)){
	prev_successful_degree = i;
	return p;
      }
      //if(p) rlibm_validate_and_fix_intervals(sp_intervals, sp_size, p);
      if(p != nullptr){
	free(p);
      }
      count++;
    }    
  }
  return nullptr;

}

int sample_compare(const void* s1, const void* s2){

  sample_info* e1 = (sample_info*) s1;
  sample_info* e2 = (sample_info*) s2;
  return e1->key > e2->key;
}

void rlibm_print_sample(sample_info* sampled_indices, size_t size){

  double prev = 0.0;
  for(size_t i = 0; i < size; i++){
    assert(sampled_indices[i].key >= prev);
    prev = sampled_indices[i].key;
    printf("counter=%lu, key=%e, sample_index=%lu\n", i, sampled_indices[i].key,
	   sampled_indices[i].index);
  }
}

void rlibm_weighted_random_sample(sample_info* sampled_indices, size_t ssize,
				  interval_data* intervals, size_t nentries){

  for(size_t i = 0; i < ssize; i++){
    sampled_indices[i].key = pow(intervals[i].u, 1./(intervals[i].w));
    sampled_indices[i].index = i;
  }
  
  qsort(sampled_indices, ssize, sizeof(sample_info), sample_compare);
  //  check_sorted (sampled_indices, ssize);

  /* select the top ssize indices from the entire interval data */
  
  for(size_t i = ssize; i < nentries; i++){

    /* if the ith element is smaller than the least element in the
       sample, then nothing to do */
    size_t j= 0;

    double interval_key = pow(intervals[i].u, 1./(intervals[i].w));
    
    if(interval_key > sampled_indices[0].key){
      /* do insertion sort of the data */
      while(interval_key > sampled_indices[j].key && j < ssize) j++;

      for(size_t t=1; t < j; t++) {
	sampled_indices[t-1] = sampled_indices[t];
      }
      sampled_indices[j-1].key = interval_key;
      sampled_indices[j-1].index = i;
    }
  }
  //  check_sorted(sampled_indices, ssize);
}


size_t rlibm_compute_violated_indices(size_t* violated_indices, interval_data* intervals, size_t nentries, polynomial* poly){

  size_t num_violated_indices = 0;
  for(size_t i = 0; i < nentries; i++){
    double y = rlibm_poly_evaluation(intervals[i].x, poly);
    if( y < intervals[i].lb || y > intervals[i].ub){
      violated_indices[num_violated_indices] = i;
      num_violated_indices++;
    }
  }
  return num_violated_indices;
}

void rlibm_evaluate_and_update_weights(size_t* violated_indices, size_t num_violated_indices,
				       interval_data* intervals, size_t nentries, size_t d){
  double w_v = 0.0;
  double w_s = 0.0;

  // this can be optimized to only change the updated weights. For now
  // using an inefficient one
  for (size_t i = 0; i < nentries; i++){
    w_s = w_s + intervals[i].w;
  }

  for(size_t i = 0; i < num_violated_indices; i++){
    w_v = w_v + intervals[violated_indices[i]].w;
  }

  //doubling the weights for a lucky iteration
  if(w_v <= 2 * w_s / (9*d -1)){
    for(size_t i = 0; i < num_violated_indices; i++){
      size_t vindex = violated_indices[i];
      intervals[vindex].w  = intervals[vindex].w * 2;
    }
  }  
}

void
rlibm_regenerate_random_values_and_reset_weights(interval_data* intervals,
						 size_t nentries){

  std::random_device rd;
  std::mt19937 generator(rd());
  std::uniform_real_distribution<double> distribution(0.0, 1.0);

  for(size_t i = 0; i < nentries; i++){
    intervals[i].u = distribution(generator);
    intervals[i].w = 1.0;
  }
}

bool check_sampled_indices(sample_info* sample, sample_info* prev_sample, size_t size){

  for(size_t i =0; i < size; i++){
    if (sample[i].index != prev_sample[i].index){
      return false;
    }
  }
  return true;
}

void rlibm_print_polyinfo(polynomial* p){

  if(p->termsize == 0){
    printf("Polynomial has no terms!\n");
    exit(0);
  }

  printf("Polynomial: y=%ax^(%d)",p->coeffs[0],p->power[0]);
  for(int j=1;j<p->termsize;j++){
    printf(" + %a x^(%d)",p->coeffs[j],p->power[j]);
  }
  printf("\n");

}

int main(int argc, char** argv){

  if(argc != 6){
    printf("Usage: %s <interval file> <violated indices><config_file> <number of violated points> <function name>\n", argv[0]);
    exit(0);
  }
  VIOLATE_THRESHOLD = atoi(argv[4]);

  set_function_process(argv);

  printf("EXIT_ON_THRESHOLD is %d\n", RLIBM_EXIT_ON_THRESHOLD);

  FILE* fpv = fopen(argv[2], "r");
  fseek(fpv, 0, SEEK_END);
  unsigned long n_v_indices = ftell(fpv);
  n_v_indices /= sizeof(double);
  printf("number of violated indices = %lu\n", n_v_indices);
  fseek(fpv, 0, SEEK_SET);

  std::unordered_set<double> sp_cases;
  for(size_t m = 0; m < n_v_indices; ++m){
      double data[1];
      size_t bytes = fread(data, sizeof(double), 1, fpv);
      sp_cases.insert(data[0]);
  }
  fclose(fpv);

  unsigned long sp_size = sp_cases.size();
  
  FILE* fp = fopen(argv[1], "r");
  assert(fp != nullptr);

  /* count the number of entries */

  fseek(fp, 0, SEEK_END);
  unsigned long nentries = ftell(fp);
  nentries = nentries/(3*sizeof(double));
  nentries -= sp_cases.size();
  printf("number of intervals = %lu\n", nentries);
  fseek(fp, 0, SEEK_SET);

  interval_data* intervals = (interval_data*) calloc(nentries, sizeof(interval_data));
  interval_data* full_intervals = (interval_data*) calloc(sp_cases.size(), sizeof(interval_data));

  std::random_device rd;
  std::mt19937 generator(rd());
  std::uniform_real_distribution<double> distribution(0.0, 1.0);

  unsigned count = 0, sp_count = 0;
  for (unsigned long i = 0; i < nentries + sp_size; i++){
    double data_entry[3];
    size_t bytes = fread(data_entry, sizeof(double), 3, fp);
    if(sp_cases.count(data_entry[0]) == 0){
      intervals[count].w = 1.0;
      intervals[count].u = distribution(generator);
      intervals[count].x = data_entry[0];
      intervals[count].lb = data_entry[1];
      intervals[count].ub = data_entry[2];
      ++count;
    }
    else{
      full_intervals[sp_count].x = data_entry[0];
      full_intervals[sp_count].lb = data_entry[1];
      full_intervals[sp_count].ub = data_entry[2];
      ++sp_count;
    }
  }
  fclose(fp); 

  auto start = high_resolution_clock::now();

  // int powers[]={0,2,4,6,8,10,12,14}; // cos
  // int powers[]={0, 1,2,3,4,5,6, 7, 8, 9}; // exp2, exp10, exp
  // int powers[] = {1,3,5,7}; //sinhx
  //int powers[] = {0,2,4}; // coshx
  //int powers[] = {1,3 ,5,7}; // sinpi
  // int powers[] = {0,2,4,6}; // cospi

  FILE* powers_file = fopen(argv[3], "r");
  assert(powers_file != nullptr);
  
  int N_RLIBM_PIECES = 1;
  int powers_size = 1;
  
  int s_ret_value = fscanf(powers_file, "%d\n", &N_RLIBM_PIECES);
  assert(s_ret_value == 1);
  
  s_ret_value = fscanf(powers_file, "%d", &powers_size);
  assert(s_ret_value == 1);
  
  int* powers = (int*) calloc(powers_size, sizeof(int));
  for(int i = 0; i < powers_size; i++){
    s_ret_value = fscanf(powers_file, "%d", &powers[i]);
    assert(s_ret_value == 1);
  }
  fclose(powers_file);
  


  /* sample size */

  size_t cd = 9 * powers_size * powers_size;
  size_t samplesize = cd;

  size_t n_violated_indices = 0, n_sp_violated_indices = 0;
  size_t *violated_indices = (size_t *) calloc(nentries, sizeof(size_t));
  size_t *sp_violated_indices = (size_t *) calloc(sp_size, sizeof(size_t));

  sample_info* sampled_indices = (sample_info*) calloc(cd, sizeof(sample_info));

  size_t prev_violated_indices = 0;
  size_t matched_violated_indices = 0;

  sample_data* sampled_intervals = (sample_data *) calloc(cd, sizeof(sample_data));

  sample_data* special_intervals = (sample_data*) calloc(sp_cases.size(), sizeof(sample_data));
  
  polynomial* p = nullptr;
  size_t total_iterations = 0;

  do{
    if(p != nullptr) free(p);
    
    n_violated_indices = 0;
    
    rlibm_weighted_random_sample(sampled_indices, cd, intervals, nentries);    
    total_iterations++;
    
    for (size_t i = 0; i < cd; i++){
      size_t iindex = sampled_indices[i].index;
      
      sampled_intervals[i].x = intervals[iindex].x;
      sampled_intervals[i].lb = intervals[iindex].lb;
      sampled_intervals[i].ub = intervals[iindex].ub;
      sampled_intervals[i].orig_lb = sampled_intervals[i].lb;
      sampled_intervals[i].orig_ub = sampled_intervals[i].ub;
      sampled_intervals[i].w = intervals[iindex].w;
      sampled_intervals[i].u = intervals[iindex].u;
      sampled_intervals[i].k = sampled_indices[i].key;
    }

    for(unsigned long k = 0; k < sp_size; ++k){
        special_intervals[k].x = full_intervals[k].x;
        special_intervals[k].lb = full_intervals[k].lb;
        special_intervals[k].ub = full_intervals[k].ub;
        special_intervals[k].orig_lb = full_intervals[k].lb;
        special_intervals[k].orig_ub = full_intervals[k].ub;
    }

    /* need to implement these functions */
    p = rlibm_generate_polynomial(sampled_intervals, samplesize, special_intervals, sp_size, powers, powers_size, MAX_TRIES);

    if(p){
      n_violated_indices = rlibm_compute_violated_indices(violated_indices, intervals, nentries, p);
      n_sp_violated_indices = rlibm_compute_violated_indices(sp_violated_indices, full_intervals, sp_size, p);
      printf("number of violated intervals: %lu, number of special violated indices: %lu, total iterations=%lu \n", n_violated_indices, n_sp_violated_indices, total_iterations);

      if(n_violated_indices+n_sp_violated_indices <= VIOLATE_THRESHOLD){
	printf("VIOLATING INPUTS BELOW THRESHOLD:\n");
	for(size_t m = 0; m < n_violated_indices; m++){
	  printf("violated_input is %a, lb is %a, ub is %a\n", intervals[violated_indices[m]].x, intervals[violated_indices[m]].lb, intervals[violated_indices[m]].ub);
	}
	for(size_t m = 0; m < n_sp_violated_indices; m++){
	  printf("sp_violated_input is %a, lb is %a, ub is %a\n", full_intervals[sp_violated_indices[m]].x, full_intervals[sp_violated_indices[m]].lb, full_intervals[sp_violated_indices[m]].ub);
	}
	rlibm_print_polyinfo(p);
	if(RLIBM_EXIT_ON_THRESHOLD){
	  break;
	}
      }
      
      rlibm_evaluate_and_update_weights(violated_indices, n_violated_indices, intervals, nentries, powers_size);

    }
    else {
      if(total_iterations > MAX_ITERATIONS){
	printf("total iterations exceeded %d, terminating the polynomial geenerator\n", MAX_ITERATIONS);
	if(p!= nullptr){
	  free(p);
	  p = nullptr;	  
	}
	break;
      }
      printf("failed to generate polynomial, resetting weights, total_iterations=%lu\n", total_iterations);
      //prev_successful_degree = 0;      
      rlibm_regenerate_random_values_and_reset_weights(intervals, nentries);
    }

    /* debugging feature to reset weights for the sample if not making progress*/
    if(n_violated_indices != 0 && (prev_violated_indices == n_violated_indices)){
      matched_violated_indices++;
      if(matched_violated_indices > SAMPLE_MATCH_THRESHOLD){
	matched_violated_indices = 0;
	n_violated_indices = 0;
	
	printf("not making progress, same number of violated indices, resetting weights, total_iterations=%lu\n", total_iterations);
	prev_successful_degree = 0;
	rlibm_regenerate_random_values_and_reset_weights(intervals, nentries);
	if(p!= nullptr) {
	  free(p);
	  p = nullptr;
	}
	continue;
      }
    }
    else{
      matched_violated_indices = 0;
      prev_violated_indices = n_violated_indices;
    }    
  } while(n_violated_indices > 0 || !p);

  if(p){
    rlibm_print_polyinfo(p);
  }
  else {
    printf("Could not generate the polynomial that satisifies all intervals, check for partial results with a few violated intervals\n");
  }

  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<minutes>(stop - start);
  printf("Total execution time: %ld\n", duration.count());

  free(p);
  free(sampled_intervals);
  free(special_intervals);
  free(sampled_indices);
  free(intervals);
  free(violated_indices);
  
  return 0;

}
