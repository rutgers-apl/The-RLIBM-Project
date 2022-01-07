#ifndef __RLIBM_COMMON_H__
#define __RLIBM_COMMON_H__

#include <stdint.h>

typedef union {
  double d;
  uint64_t x;
} double_x;

typedef union {
  float f;
  uint32_t x;
} float_x;


#endif
