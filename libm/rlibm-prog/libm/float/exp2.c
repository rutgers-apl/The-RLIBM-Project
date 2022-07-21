#include <math.h>
#include "float_lib.h"
#include "exp2.h"

float rlibm_fast_exp2(float x) {
    float_x fx;
    fx.f = x;
    
    // Take care of special cases
    if (0x43000000 <= fx.x && fx.x <= 0xb338aa3b) {
        if (fx.x <= 0x7F800000) return 1.0/0.0;
        if (fx.x < 0x80000000) return 0.0/0.0;
        return 1.0;
    }
    
    if (fx.x <= 0x33b8aa3a) {
        return 1.0;
    }
    
    if (fx.x >= 0xc3160000) {
        if (fx.x <= 0xFF800000) return 0.0;
        return 0.0/0.0;
    }
    
    // Perform range reduction
    double xp = x * 64;
    int N = (int)xp;
    int N2 = N % 64;
    if (N2 < 0) N2 += 64;
    int N1 = N - N2;
    
    int M = N1 / 64;
    int J = N2;
    double R = x - N * 0.015625;

    double y = 0.0;
    
    y = 1.3529558545669370488162552845778918708674609661102294921875000000000000e-03;
    y *= R;
    y += 9.6181130451127132274802278288916568271815776824951171875000000000000000e-03;
    y *= R;
    y += 5.5504104858622581308846832826020545326173305511474609375000000000000000e-02;
    y *= R;
    y += 2.4022650697767972127749658284301403909921646118164062500000000000000000e-01;
    y *= R;
    y += 6.9314718055997448509231162461219355463981628417968750000000000000000000e-01;
    y *= R;
    y += 9.9999999999999988897769753748434595763683319091796875000000000000000000e-01;

#if 0    
    
    double_x dX;        // TODO: Are these 2 lines redundant?
    dX.d = R;

    unsigned long Rbits = *(unsigned long*)&R;
    double y = 1.3529558545669370488162552845778918708674609661102294921875000000000000e-03;

    if(Rbits & 0xBFFFFFFFE0000000 == Rbits){  // satisfied by 1383939 out of 302951375 intervals.
        y *= R;
        y += 9.6181130451127132274802278288916568271815776824951171875000000000000000e-03;
        y *= R;
        y += 5.5504104858622581308846832826020545326173305511474609375000000000000000e-02;
        y *= R;
        y += 2.4022650697767972127749658284301403909921646118164062500000000000000000e-01;
        y *= R;
        y += 6.9314718055997448509231162461219355463981628417968750000000000000000000e-01;
        y *= R;
        y += 9.9999999999999988897769753748434595763683319091796875000000000000000000e-01;
    }
    else{
        y = 5.6635361177615944905383571494894567877054214477539062500000000000000000e-02;
        y *= R;
        y += 2.4022841172163697520680614161392441019415855407714843750000000000000000e-01;
        y *= R;
        y += 6.9314697903493327491020181696512736380100250244140625000000000000000000e-01;
        y *= R;
        y += 1.0000000000298172597723578292061574757099151611328125000000000000000000e+00;
    }
#endif    
    // Perform output compensation
    return y * ldexp(exp2JBy64[J], M);
}