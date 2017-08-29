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
  float fa = 25235.f;
  float fb = -0.53123f;
  float32 a = fa;
  float32 b = fb;
  std::cout<<float(a)<<'\n';
  if(a == b) {std::cout<<"SUCCESS\n";}
  auto c = a+b;
  std::printf("val: %.64f\n",float(c));
  std::printf("val: %.64f\n",float(fa+fb));
  return 0;
}
