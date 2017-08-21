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
  double a = 0.5;
  float64 b(a);
  float16 z(b);
  std::printf("val: %.64f\n",double(b));
  std::printf("val: %.64f\n",double(z));
  std::cout<<std::numeric_limits<float32>::is_specialized<<'\n';
  std::printf("val: %.64f\n",float(std::numeric_limits<float64>::max_exponent10));
  std::printf("val: %.64f\n",float(std::numeric_limits<double>::max_exponent10));
  std::cout<<std::log2(std::numeric_limits<float>::max())<<'\n';
  int bar;
  std::frexp(std::numeric_limits<float>::max(),&bar);
  std::cout<<(1<<7)/std::log2(10)<<'\n';
  auto foo = bit_cast<uint64_t>(std::numeric_limits<double>::max());
  std::cout<<std::bitset<64>(bar)<<'\n';
  return 0;
}
//std::numeric_limits<T>::is_specialized
