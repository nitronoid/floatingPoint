#include <cstdlib>
#include <iostream>
#include <cmath>
#include <bitset>
#include "dynamicfloat.h"

using float8 = dynamicFloat<3,4>;
using float16 = dynamicFloat<10,5>;
using float32 = dynamicFloat<23,8>;
using float64 = dynamicFloat<52,11>;

int main()
{
  float a = 10000000000000000000;
  float32 b(a);
  std::printf("val: %.64f\n",float(a));
  std::printf("val: %.64f\n",float(std::numeric_limits<float>::denorm_min()));
  return 0;
}
//std::numeric_limits<float32>::epsilon().getBits() (uint64_t(1)<<(N-1))-1 - ((uint64_t(1)<<MANTISSA)-1)
