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
  float a = 0.1;
  //float8 x(a);
  float32 z(a);
  std::printf("val: %.64f\n",a);
  //std::printf("val: %.64f\n",double(x));
  std::printf("val: %.64f\n",float(z));
  return 0;
}
