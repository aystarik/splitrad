#define NEED_TEST_DATA
#include "fft.h"
int main()
{
  float cshift, rshift;
  test(cshift, rshift);
  printf("%f,%f\n", static_cast<double>(cshift), static_cast<double>(rshift));
  return 0;
}
