#include "float-rno/float_rno_lib.h"
#include "log10.h"
#include "math.h"

#define LOG102HIGH 0.30102999566398114250631579125183634459972381591796875
#define LOG102LOW  5.27074231034726570126349709198449199648263806413338306011695522101945243775844573974609375e-17

double rlibm_rno_log10(float x) {
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

  switch (fix.x) {
  case 0x3f800000 : return 0.0;
  case 0x41200000 : return 1.0;
  case 0x42c80000 : return 2.0;
  case 0x447a0000 : return 3.0;
  case 0x461c4000 : return 4.0;
  case 0x47c35000 : return 5.0;
  case 0x49742400 : return 6.0;
  case 0x4b189680 : return 7.0;
  case 0x4cbebc20 : return 8.0;
  case 0x4e6e6b28 : return 9.0;
  case 0x501502f9 : return 10.0;
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
  
  // Find the index of polynomial coefficients
  double_x dX;
  dX.d = f;
  
  double y;

  switch (dX.x) {
    case 0x3f44f22d0e560419:
      y = 2.775213025069361958037383875108616848592646420001983642578125e-04;
      break;
    case 0x3f61fddb0d3224f3:
      y = 9.52769693548433972696276583747021504677832126617431640625e-04;
      break;
    case 0x3f68ff099fc267f0:
      y = 1.32314440169836928552771215805705651291646063327789306640625e-03;
      break;
    case 0x3f6bde34a2b10bf6:
      y = 1.4748992865482070391269786568955169059336185455322265625e-03;
      break;
    case 0x3f6fbf5f5f5f5f5f:
      y = 1.679826337715653768178913907149762962944805622100830078125e-03;
      break;
    case 0x3f72af84a062b2e5:
      y = 1.97671055057299822899086194638584856875240802764892578125e-03;
      break;
    case 0x3f74212f684bda13:
      y = 2.1290956120136798716824255706114854547195136547088623046875e-03;
      break;
    default :
      y = 9.6039176808723147882318471602047793567180633544921875e-02;
      y *= f;
      y += -1.087281288826748848475034492366830818355083465576171875e-01;
      y *= f;
      y += 1.447656742943796148725965622361400164663791656494140625e-01;
      y *= f;
      y += -2.171472428171018209663856168845086358487606048583984375e-01;
      y *= f;
      y += 4.34294481904594020793552999748499132692813873291015625e-01;
      y *= f;
  }
  
  y += m * LOG102LOW;
  y += log10_lut[FIndex];
  y += m * LOG102HIGH;
  
  return y;
}
