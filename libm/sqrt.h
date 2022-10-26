#if defined(__SSE2__)

#include <emmintrin.h>

static double Sqrt(double x)
{
    __m128d x_vec = _mm_set_pd(x, x);
    __m128d sqrt_x = _mm_sqrt_sd(x_vec, x_vec);
    return _mm_cvtsd_f64(sqrt_x);
}

#else

#include <math.h>

static double Sqrt(double x)
{
  return sqrt(x);
}

#endif
