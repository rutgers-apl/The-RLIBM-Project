#include <math.h>
#include <mpfr.h>
#include <stdio.h>
#include<stdbool.h>

#include "float34rno.h"


mpfr_t mval, r;

int main(int argc, char** argv) {
  mpfr_inits(mval, r, (mpfr_ptr) 0);
  mpfr_set_prec(mval, 1000);
  /* E Min = (-126) - 25 + 1 = -150 */
  /* E Max = 127 + 1 = 128 */
  mpfr_set_emin(-150);
  mpfr_set_emax(128);
  float_x x;
  double_x y;
  FILE *fp;
  fp = fopen(argv[1], "a+");
  for (unsigned long count = 0; count < 0x100000000; count++) {
    if(count % 0x100000 == 0){
      printf("Completed count = %lx \n", count);
    }
    if ((count & 0x7FFFFFFF) > 0x7F800000) {
      y.d =  0.0f/0.0f;
      fwrite(&y.d, sizeof(double), 1, fp);
      continue;
    }
    x.x = count;
    /* Set mval to be the floating point value of input x.*/
    int sticky = mpfr_set_d(mval, x.f, MPFR_RNDZ);
    /* Sticky should be 0 up to this point. (add assert) */
    /* Compute sin(x) and subnormalize. */
    sticky = mpfr_sin(mval, mval, MPFR_RNDZ); 
    sticky = mpfr_subnormalize(mval, sticky, MPFR_RNDZ);
    y.d = getFloat34RNO(mval, r, sticky);
    fwrite(&y.d, sizeof(double), 1, fp);
  }
  fclose(fp);
  mpfr_clears(mval, r, (mpfr_ptr) 0);
  return 0;
}

