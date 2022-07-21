#include "sinf.h"

int main(int argc, char** argv){
  float inp = 1.572;

  float result = rlibm_fast_sin(inp);
  printf("the result is %a\n", result);

  return 0;

}
