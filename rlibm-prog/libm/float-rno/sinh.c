#include "float-rno/float_rno_lib.h"
#include "float_lib.h"
#include "float-rno/sinhcosh.h"

#define CONST64BYLN2 92.332482616893656768297660164535045623779296875
#define LN2BY64 0.01083042469624914509729318723429969395510852336883544921875

double rlibm_rno_sinh(float x) {
  float_x fx;
  fx.f = x;
  unsigned long sign = (fx.x & 0x80000000) == 0 ? 0x0 : 0x8000000000000000;
  fx.x &= 0x7FFFFFFF;

  if (fx.x == 0) return 0.0;
  
  // Take care of special cases
  if (fx.x <= 971544424) {
    double_x dX;
    dX.d = (double)fx.f;
    long exp = (dX.x & 0x7FF0000000000000UL) >> 52UL;
    exp -= 1023L;
    long mantissaCount = exp + 149L;
    if (mantissaCount > 23) mantissaCount = 23;
    mantissaCount += 2L;
    unsigned long shiftAmount = (52L - mantissaCount);
    unsigned long sticky = 1UL << shiftAmount;
    dX.x |= sticky;
    dX.x |= sign;
    return dX.d;
  }
  
  if (fx.x >= 1119016189) {
    if (fx.x > 0x7F800000) return 0.0/0.0;
    if (fx.x == 0x7F800000) {
      if (x > 0.0f) return 1.0 / 0.0;
      else return -1.0 / 0.0;
    }

    if (x > 0.0f) return 3.40282361850336062550457001444955389952e+38;
    else return -3.40282361850336062550457001444955389952e+38;
  }
  
  // Perform range reduction
  double xp = fx.f * CONST64BYLN2;
  int N = (int)xp;
  int N2 = N % 64;
  if (N2 < 0) N2 += 64;
  int N1 = N - N2;
  int I = N1 / 64;
  double R = fx.f - N * LN2BY64;
  double R2 = R * R;
  
  double sinhHigh = sinhKLn2[I];
  double coshHigh = coshKLn2[I];
  double sinhMid = sinhKLn2By64[N2];
  double coshMid = coshKLn2By64[N2];
  
  double sinhHM = sinhHigh * coshMid + coshHigh * sinhMid;
  double coshHM = sinhHigh * sinhMid + coshHigh * coshMid;
  
  // Compute sinh component
  double sinhL;
  if (R == 8.74292667205892222448415651570030604489147663116455078125e-03) {
    sinhL = 8.7430380555734432679315659697749651968479156494140625e-03;
  } else if (R == 1.03883635115007422200505970977246761322021484375e-02) {
    sinhL = 1.0388550361244829056683869339394732378423213958740234375e-02;
  } else if (R == 1.0487717054569856145462836138904094696044921875e-02) {
    sinhL = 1.0487909316821579508438588845820049755275249481201171875e-02;
  } else if (R == 1.0769002139568328857421875e-02) {
    sinhL = 1.076921029016375715159359316430709441192448139190673828125e-02;
  } else {
    sinhL = 8.3328876580867176915301541839653509669005870819091796875e-03;
    sinhL *= R2;
    sinhL += 1.666666666736031088280611811569542624056339263916015625e-01;
    sinhL *= R2;
    sinhL += 9.999999999999997779553950749686919152736663818359375e-01;
    sinhL *= R;
  }
  
  double coshL;
  if (R == 8.3387088168791478892671875655651092529296875e-03) {
    coshL = 1.000034736701873594455491911503486335277557373046875;
  } else {
    coshL = 4.1667466122972861286566370608852594159543514251708984375e-02;
    coshL *= R2;
    coshL += 4.99999999948495255086555744128418155014514923095703125e-01;
    coshL *= R2;
    coshL += 1.0000000000000002220446049250313080847263336181640625;
  }
  
  // Perform output compensation
  double_x dX;
  dX.d = sinhHM * coshL + coshHM * sinhL;
  dX.x |= sign;
  return dX.d;
}
