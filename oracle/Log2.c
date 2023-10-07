#include "OracleGenerate.h"

double ComputeOracleResult(float x) {
    if (x == 1.0 / 0.0) return 1.0 / 0.0;
    if (x == 0) return -1.0 / 0.0;
    
    int sticky = 0;
    // Set float value to mpfr. This should be exact
    int status = mpfr_set_d(mval, (double)x, MPFR_RNDN);
    if (status != 0) {
        printf("Something went wrong when setting float to mpfr\n");
        exit(0);
    }
    
    // Call the elementary function. Make it round to zero. Then we can tell
    // what the sticky bit is from the return value.
    status = mpfr_log2(mval, mval, MPFR_RNDZ);
    if (status != 0) sticky |= 0x1;
    
    return FromMPFRToFloat34Ro(mval, sticky);
}

int main(int argc, char** argv) {
    mpfr_init2(mval, 200);

    if (argc != 2) {
        printf("UsageL %s <output file>\n", argv[0]);
        exit(0);
    }
    
    RunTest(argv[1]);
    mpfr_clear(mval);
    return 0;
}
