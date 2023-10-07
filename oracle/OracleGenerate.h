#include "stdio.h"
#include "stdlib.h"
#include "RoundToOdd.h"

mpfr_t mval;

double ComputeOracleResult(float x);

void RunTestHelper(char* FileName, unsigned long start, unsigned long end) {
  unsigned long count = 0;
  
  float x;
  floatX xbase;
  
  FILE* f = fopen(FileName, "w");
  
  for (count = start; count < end; count++) {
    xbase.x = count;
    x = xbase.f;
    
    double res = ComputeOracleResult(x);

    if(count % 0x100000 == 0){
      printf("Completed count = %lx \n", count);
    }        
    fwrite(&res, sizeof(double), 1, f);
  }
  
  fclose(f);

}

void RunTest(char* FileName) {
  RunTestHelper(FileName, 0x0, 0x100000000);
}
