#include<stdio.h>
#include<stdlib.h>
#include "sinf.h"

int main(int argc, char** argv){

  FLT x;

  unsigned long count = 0;
  unsigned long wrong_count = 0;

  FILE* fp = fopen(argv[1], "r");

  for (count = 0; count < 0x100000000; count++){

    if(count % 10000000 == 0){
      printf("Completed count = %lx \n", count);
    }
    x.i = count;
    FLT result;
    result.f = rlibm_fast_sin(x.f);
    FLT oracle;
    fread(&oracle.f, sizeof(float), 1, fp);
    if(oracle.i != result.i){
      printf("mismatch\n");
      wrong_count++;
    }
    
  }
  fclose(fp);
  if(wrong_count == 0){
    printf("Correct results for all inputs\n");
  }
  else{
    printf("Wrong results for %lu inputs\n", wrong_count);
  }
  return 0;
}
