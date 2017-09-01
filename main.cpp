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
  double fa = 5460.4365;
  double fb = 206456.44455;
  float32 a = fa;
  float32 b = fb;
  std::cout<<float(a)<<'\n';
  if(a == b) {std::cout<<"SUCCESS\n";}
  auto c = a*b;
  std::printf("val: %.64f\n",float(c));
  std::printf("val: %.64f\n",float(fa*fb));
  std::cout<<std::bitset<32>(bit_cast<uint32_t>(float(c)))<<'\n';
  std::cout<<std::bitset<32>(bit_cast<uint32_t>(float(fa*fb)))<<'\n';
  return 0;
}
