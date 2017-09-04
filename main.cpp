#include <cstdlib>
#include <iostream>
#include <cmath>
#include <bitset>
#include "dynamicfloat.h"

using float8 = dynamicFloat<3,4>;
using float16 = dynamicFloat<10,5>;
using float32 = dynamicFloat<23,8>;
using float64 = dynamicFloat<52,11>;

#define OPERATION /

int main()
{
  using test = double;
  using testW = uint64_t;
  using testDf = float64;
  constexpr int width = 64;

  test fa = 486759365936596563;
  test fb = 34.534;
  testDf a = fa;
  testDf b = fb;

  if(a == b) {std::cout<<"SUCCESS\n";}
  auto c = a OPERATION b;
  std::printf("val: %.64f\n",test(c));
  std::printf("val: %.64f\n",test(fa OPERATION fb));
  std::cout<<std::bitset<width>(bit_cast<testW>(test(c)))<<'\n';
  std::cout<<std::bitset<width>(bit_cast<testW>(test(fa OPERATION fb)))<<'\n';
  std::cout<<sizeof(long double)<<'\n';
  return 0;
}
