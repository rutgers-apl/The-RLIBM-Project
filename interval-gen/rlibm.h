#pragma once

#include "mpfr.h"

using namespace std;

typedef union {
    float f;
    unsigned x;
} float_x;

typedef union {
    double d;
    unsigned long long int x;
} double_x;

struct IntData {
    double lb;
    double ub;
};
