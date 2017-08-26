#ifndef DYNAMICFLOAT_H
#define DYNAMICFLOAT_H

#include <cstdint>
#include <cmath>
#include <iostream>
#include <limits>
#include <cstring>

template <typename DEST, typename SOURCE>
inline DEST bit_cast(const SOURCE& source)
{
  static_assert(sizeof(DEST) == sizeof(SOURCE),"Dest and source must have same size!");
  DEST dest;
  std::memcpy(&dest, &source, sizeof(dest));
  return dest;
}

// Selects the correct size type based on passed bit num
template<int N>
using uintN_t =
  typename std::conditional_t< N <= 8, uint8_t,
    typename std::conditional_t< N <= 16, uint16_t,
      typename std::conditional_t< N <= 32, uint32_t,
        uint64_t
    >
  >
>;

template <uint64_t SIGNIFICAND, uint64_t EXPONENT>
class dynamicFloat
{
  //Friend the numeric_limits specialisation
  friend struct std::numeric_limits<dynamicFloat<SIGNIFICAND, EXPONENT>>;
  //Friend all specialisations of dynamic float for easier conversion
  template<uint64_t,uint64_t > friend class dynamicFloat;
  //Declare total number of bits for easy reference
  static constexpr int typeWidth = SIGNIFICAND + EXPONENT + 1;
  //Define the smallest unsigned int that can represent this
  using uint_n = uintN_t<typeWidth>;

public:

  dynamicFloat() = default;
  explicit dynamicFloat(float _infloat) : dynamicFloat(bit_cast<dynamicFloat<23,8>>(_infloat)){}
  explicit dynamicFloat(double _indouble) : dynamicFloat(bit_cast<dynamicFloat<52,11>>(_indouble)){}

  template<uint64_t M, uint64_t E>
  dynamicFloat(const dynamicFloat<M,E> &_inFloat)
  {
    using iFloat_t = dynamicFloat<M,E>;
    // sign bit should remain the same
    IEEE.sign = _inFloat.IEEE.sign;

    // calculate the exponent value without bias, single has 127 bias
    int64_t unbiasedExponent = _inFloat.IEEE.exponent-bias(iFloat_t::exponentLength());

    // find the difference in size between the mantissas
    int64_t padSize = padLen<dynamicFloat<M,E>>();

    if(unbiasedExponent < minExponent())
    {
      IEEE.significand = IEEE.exponent = 0;
    }

    // if NaN or Infinity, stay NaN or Infinity
    else if (_inFloat.IEEE.exponent == exponentMaskOnes(iFloat_t::exponentLength()))
    {
      // NaN or INF
      IEEE.significand = _inFloat.IEEE.significand;
      IEEE.exponent = exponentMaskOnes(EXPONENT);
    }

    else if (unbiasedExponent < (1 - bias(EXPONENT)))
    {
      // this maps to a subnormal number
      IEEE.exponent = 0;
      int64_t shiftVal =  1 - (unbiasedExponent + bias(EXPONENT));
      IEEE.significand = (base(SIGNIFICAND) >> shiftVal) + ( _inFloat.IEEE.significand >> (padSize + shiftVal));
    }

    // if the exponent is zero, map to zero
    else if (!_inFloat.IEEE.exponent)
    {
      if(!_inFloat.IEEE.significand)
        IEEE.significand = IEEE.exponent = 0;
      else
      {
        uint64_t e = 0;
        uint64_t m = _inFloat.IEEE.significand << 1;

        for(;!(m & base(iFloat_t::significandLength())); ++e, m <<= 1);

        IEEE.exponent = bias(EXPONENT) - bias(iFloat_t::exponentLength()) - e;
        IEEE.significand = (padSize < 0) ?
              ((m & (base(iFloat_t::significandLength()) - 1)) << -padSize )
            : ((m & (base(iFloat_t::significandLength()) - 1)) >> padSize);
      }
    }

    else if (unbiasedExponent > bias(EXPONENT))
    {
      // too large to be represented, map this value to infinity
      IEEE.significand = 0;
      IEEE.exponent = exponentMaskOnes(EXPONENT);
    }

    else
    {
      // normal numbers
      IEEE.exponent = unbiasedExponent + bias(EXPONENT); // add the new bias
      uint64_t signicandCopy = _inFloat.IEEE.significand;
      IEEE.significand = (padSize < 0) ? (signicandCopy << -padSize) : (signicandCopy >> padSize);
    }
  }

  operator float() const
  {
    return bit_cast<float>(dynamicFloat<23,8>(*this));
  }

  operator double() const
  {
    return bit_cast<double>(dynamicFloat<52,11>(*this));
  }

private:
  // union to store this classes data
  union
  {
    uint_n m_bits;                // All bits
    struct
    {
      uint_n significand : SIGNIFICAND;	// mantissa
      uint_n exponent : EXPONENT;	// exponent
      uint_n sign : 1;            // sign always 1 bit
    } IEEE;
  };

  static constexpr int64_t bias(uint64_t _exponent) noexcept { return (uint64_t{1} << (_exponent-1)) - 1; }

  static constexpr int64_t minExponent() noexcept { return -(SIGNIFICAND - 1 + bias(EXPONENT)); }

  static constexpr uint64_t exponentMaskOnes(uint64_t _exponent) noexcept { return (uint64_t{1} << _exponent) - 1; }

  static constexpr uint64_t base(uint64_t _mantissa) noexcept { return uint64_t{1} << _mantissa; }

  static constexpr uint64_t significandLength() noexcept { return SIGNIFICAND; }

  static constexpr uint64_t exponentLength() noexcept { return EXPONENT; }

  template<typename T>
  static constexpr int64_t padLen() noexcept { return T::significandLength() - SIGNIFICAND; }
};

namespace std
{
template <uint64_t SIGNIFICAND, uint64_t EXPONENT>
struct numeric_limits<dynamicFloat<SIGNIFICAND, EXPONENT>>
{
  using dfloat = dynamicFloat<SIGNIFICAND, EXPONENT>;
  static constexpr int typeWidthMinusOne = SIGNIFICAND + EXPONENT;
  static constexpr int typeWidth   = typeWidthMinusOne + 1;
  using uint_n = uintN_t<typeWidth>;
  static constexpr uint_n maskOnes = ~uint_n{0};
  static constexpr uint_n maskOnesNoSign = maskOnes & ~(uint64_t{1} << typeWidthMinusOne); // 0-111.... , 0 followed by all ones
  static constexpr uint64_t mantAllOne(int64_t _offset) { return ((uint64_t{1} << (SIGNIFICAND + _offset)) - 1 ); }
  static constexpr dfloat floatFromBits(uint_n _bits)
  {
    //Default constructor will avoid non-constexpr call to bit_cast
    dfloat bitFloat{};
    //Set the bits directly with friend class relationship
    bitFloat.m_bits = _bits;
    return bitFloat;
  }
  static constexpr double log_base2_10 = 3.32192809489;
  static constexpr double log_base10_radix = 0.30102999566;

public:

  // General -- meaningful for all specializations.
  static constexpr bool is_specialized = true;
  static constexpr dfloat min() { return floatFromBits(uint_n{mantAllOne(0) + 1}); }
  static constexpr dfloat max() { return floatFromBits(uint_n{maskOnesNoSign - (mantAllOne(0) + 1)}); }
  static constexpr dfloat lowest() { return floatFromBits(uint_n{maskOnes - (mantAllOne(0) + 1)}); }
  static constexpr int  radix        = 2;
  static constexpr int  digits       = SIGNIFICAND + 1;
  static constexpr int  digits10     = SIGNIFICAND * log_base10_radix;
  static constexpr int  max_digits10 = digits * log_base10_radix + 2;
  static constexpr bool is_signed		 = true;
  static constexpr bool is_integer	 = false;
  static constexpr bool is_exact		 = false;
  static constexpr bool traps        = false;
  static constexpr bool is_modulo		 = false;
  static constexpr bool is_bounded	 = true;

  // Floating point specific.

  static constexpr dfloat epsilon()     noexcept { return dfloat(1.0f / (1 << SIGNIFICAND)); }
  static constexpr dfloat round_error() noexcept { return dfloat(0.5); }

  static constexpr int max_exponent   = 1 << (EXPONENT - 1);
  static constexpr int min_exponent   = 3 - max_exponent;

  static constexpr int min_exponent10 = min_exponent / log_base2_10;
  static constexpr int max_exponent10 = max_exponent / log_base2_10; //divide by log2(10)

  static constexpr bool has_infinity        = true;
  static constexpr bool has_quiet_NaN       = true;
  static constexpr bool has_signaling_NaN		= true;
  static constexpr bool is_iec559           = false;
  static constexpr bool has_denorm          = denorm_present;
  static constexpr bool tinyness_before     = false;
  static constexpr float_round_style round_style = round_to_nearest;

  static constexpr dfloat denorm_min()    noexcept { return floatFromBits(uint_n{1}); }
  static constexpr dfloat infinity()      noexcept { return floatFromBits(uint_n{maskOnesNoSign - mantAllOne(0)}); }
  static constexpr dfloat quiet_NaN()     noexcept { return floatFromBits(uint_n{maskOnesNoSign - mantAllOne(-1)}); }
  static constexpr dfloat signaling_NaN() noexcept { return floatFromBits(uint_n{maskOnesNoSign - mantAllOne(-1) - mantAllOne(-2)}); }
};
}

#endif // DYNAMICFLOAT_H
