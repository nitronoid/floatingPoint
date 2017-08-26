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
  double a = 321.435;
  float64 b(a);
  float16 z(b);
  std::printf("val: %.64f\n",double(b));
  std::printf("val: %.64f\n",double(z));
  std::cout<<float(std::numeric_limits<float32>::lowest())<<'\n';
  uint32_t foo = bit_cast<uint32_t>(std::numeric_limits<float32>::lowest());
  std::cout<<std::bitset<32>(foo)<<'\n';
  std::cout<<std::numeric_limits<float32>::has_infinity<<'\n';
  return 0;
}
//std::numeric_limits<T>::is_specialized
