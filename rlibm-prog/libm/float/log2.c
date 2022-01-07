#include <math.h>
#include "float_lib.h"
#include "log2.h"

float rlibm_fast_log2(float x) {
    float_x fix, fit;
    fix.f = x;
    int m = 0;
    
    if (fix.x < 0x800000 || fix.x >= 0x7F800000) {
        if ((fix.x & 0x7FFFFFFF) == 0) { // log(+/-0) = -infty
            fix.x = 0xFF800000;
            return fix.f;
        }
        
        if (fix.x > 0x7FFFFFFF) { // Log(-val) = NaN
            return (x - x) / 0;
            
        }
        
        if (fix.x >= 0x7F800000) {
            return x + x;
        }
        
        fix.f *= 8.388608e+06;
        m -= 23;
    }
    
    m += fix.x >> 23;
    m -= 127;
    fix.x &= 0x007FFFFF;
    fix.x |= 0x3F800000;
    
    fit.x = fix.x & 0x007F0000;
    int FIndex = fit.x >> 16;
    fit.x |= 0x3F800000;
    
    double f = fix.f - fit.f;
    f *= log_oneByF[FIndex];


    double y = 3.8662151480904477507394290114461909979581832885742187500000000000000000e-01;
    y *= f;
    y += -3.6229559281288897798489756496564950793981552124023437500000000000000000e-01;
    y *= f;
    y += 4.8090713122852124516981575652607716619968414306640625000000000000000000e-01 ;
    y *= f;
    y += -7.2134753829891251619699232833227142691612243652343750000000000000000000e-01 ;
    y *= f;
    y += 1.4426950408983432172504990376182831823825836181640625000000000000000000e+00 ;
    y *= f;

    
    return y + log2_lut[FIndex] + m;
}

