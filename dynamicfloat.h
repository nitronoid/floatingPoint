#ifndef DYNAMICFLOAT_H
#define DYNAMICFLOAT_H

#include <cstdint>
#include <cmath>
#include <iostream>
#include <limits>
#include <cstring>

template <typename DEST, typename SOURCE>
inline DEST bit_cast(const SOURCE& source) {
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

template <uint64_t MANTISSA, uint64_t EXPONENT>
class dynamicFloat
{
  //Friend the numeric_limits specialisation
  template<class> friend class std::numeric_limits;
  //Friend all specialisations of dynamic float for easier conversion
  template<uint64_t,uint64_t > friend class dynamicFloat;
  //Declare total number of bits for easy reference
  static constexpr int DIGITS = MANTISSA + EXPONENT + 1;
  //Define the smallest unsigned int that can represent this
  using uint_n = uintN_t<DIGITS>;

public:

  dynamicFloat() = default;
  dynamicFloat(float _infloat) : dynamicFloat(bit_cast<dynamicFloat<23,8>>(_infloat)){}
  dynamicFloat(double _indouble) : dynamicFloat(bit_cast<dynamicFloat<52,11>>(_indouble)){}

  template<uint64_t M, uint64_t E>
  dynamicFloat(const dynamicFloat<M,E> &_infloat)
  {
    using iFloat_t = dynamicFloat<M,E>;
    // sign bit should remain the same
    IEEE.sign = _infloat.IEEE.sign;

    // calculate the exponent value without bias, single has 127 bias
    int64_t unbiasedExponent = _infloat.IEEE.exponent-bias(iFloat_t::expLen());

    // find the difference in size between the mantissas
    int64_t padSize = padLen<dynamicFloat<M,E>>();

    if(unbiasedExponent < minExponent())
    {
      IEEE.mantissa = IEEE.exponent = 0;
    }

    // if NaN or Infinity, stay NaN or Infinity
    else if (_infloat.IEEE.exponent == maxExp(iFloat_t::expLen()))
    {
      // NaN or INF
      IEEE.mantissa = _infloat.IEEE.mantissa;
      IEEE.exponent = maxExp(EXPONENT);
    }

    else if (unbiasedExponent<1-bias(EXPONENT))
    {
      // this maps to a subnormal number
      IEEE.exponent = 0;
      int64_t shiftVal =  1-(unbiasedExponent + bias(EXPONENT));
      IEEE.mantissa = (base(MANTISSA) >> shiftVal) + ( _infloat.IEEE.mantissa >> (padSize + shiftVal));
    }

    // if the exponent is zero, map to zero
    else if (!_infloat.IEEE.exponent)
    {
      if(!_infloat.IEEE.mantissa)
        IEEE.mantissa = IEEE.exponent = 0;
      else
      {
        uint64_t e = 0;
        uint64_t m = _infloat.IEEE.mantissa<<1;

        for(;!(m & base(iFloat_t::mantLen())); ++e, m <<= 1);

        IEEE.exponent = bias(EXPONENT) - bias(iFloat_t::expLen()) - e;
        IEEE.mantissa = (padSize < 0) ? ((m & (base(iFloat_t::mantLen())-1)) << -padSize ) : ((m & (base(iFloat_t::mantLen())-1)) >> padSize);
      }
    }

    else if (unbiasedExponent > bias(EXPONENT))
    {
      // too large to be represented, map this value to infinity
      IEEE.mantissa = 0;
      IEEE.exponent = maxExp(EXPONENT);
    }

    else
    {
      // normal numbers
      IEEE.exponent = unbiasedExponent + bias(EXPONENT); // add the new bias
      uint64_t foo = _infloat.IEEE.mantissa;
      IEEE.mantissa = (padSize < 0) ? (foo << -padSize) : (foo >> padSize);
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
      uint_n mantissa : MANTISSA;	// mantissa
      uint_n exponent : EXPONENT;	// exponent
      uint_n sign : 1;            // sign always 1 bit
    } IEEE;
  };

  static constexpr int64_t bias(uint64_t _exponent) noexcept { return (uint64_t(1) << (_exponent-1)) - 1; }

  static constexpr int64_t minExponent() noexcept { return -(MANTISSA - 1 + bias(EXPONENT)); }

  static constexpr uint64_t maxExp(uint64_t _exponent) noexcept { return (uint64_t(1) << _exponent) - 1; }

  static constexpr uint64_t base(uint64_t _mantissa) noexcept { return uint64_t(1) << _mantissa; }

  static constexpr uint64_t mantLen() noexcept { return MANTISSA; }

  static constexpr uint64_t expLen() noexcept { return EXPONENT; }

  template<typename T>
  static constexpr int64_t padLen() noexcept { return T::mantLen() - MANTISSA; }
};

namespace std
{
template <int MANTISSA, int EXPONENT>
class numeric_limits<dynamicFloat<MANTISSA, EXPONENT>>
{
  using dfloat = dynamicFloat<MANTISSA, EXPONENT>;
  static constexpr int Nmo = MANTISSA + EXPONENT; // N minus one
  static constexpr int N   = Nmo + 1;
  using uint_n = uintN_t<N>;
  static constexpr uint64_t allOne = (uint64_t(1)<<Nmo)-1; // 0-111.... , 0 followed by all ones
  static constexpr uint64_t mantAllOne(int64_t _offset = 0) { return ((uint64_t(1)<<(MANTISSA + _offset)) - 1 ); }

  static constexpr uint64_t ilog2(int _n) { return ( (_n<2) ? 1 : 1+ilog2(_n/2)); }
  static constexpr uint64_t ilog10(int _n) { return ilog2(_n) / ilog2(10); }
  static constexpr dfloat floatFromBits(uint_n _bits)
  {
    //Default constructor will avoid non-constexpr call to bit_cast
    dfloat minfloat{};
    //Set the bits directly with friend class relationship
    minfloat.m_bits = _bits;
    return minfloat;
  }

public:

  // General -- meaningful for all specializations.
  static constexpr bool is_specialized = true;
  static constexpr dfloat min() { return floatFromBits(uint_n( mantAllOne() + 1 )); }
  static constexpr dfloat max() { return floatFromBits(uint_n( allOne - (mantAllOne() + 1) )); }
  static constexpr int  radix       = 2;
  static constexpr int  digits      = MANTISSA + 1;
  static constexpr int  digits10    = MANTISSA * 0.30102999566;  // *log10(radix)
  static constexpr bool is_signed		= true;
  static constexpr bool is_integer	= false;
  static constexpr bool is_exact		= false;
  static constexpr bool traps       = false;
  static constexpr bool is_modulo		= false;
  static constexpr bool is_bounded	= true;

  // Floating point specific.

  static constexpr dfloat epsilon() noexcept { return dfloat(std::pow(2,-MANTISSA)); } //Inefficient!
  static constexpr dfloat round_error() noexcept { return dfloat(0.5); }
  //static const int min_exponent10     = std::log10(double(min()));
  //static const int max_exponent10     = ilog10(max());
  static constexpr int min_exponent   = 1 - ((1 << (EXPONENT-1))-2);
  static constexpr int max_exponent   = 1<<(EXPONENT-1);

  static constexpr bool has_infinity        = true;
  static constexpr bool has_quiet_NaN       = true;
  static constexpr bool has_signaling_NaN		= true;
  static constexpr bool is_iec559           = false;
  static constexpr bool has_denorm          = denorm_present;
  static constexpr bool tinyness_before     = false;
  static constexpr float_round_style round_style = round_to_nearest;

  static constexpr dfloat denorm_min()    noexcept { return floatFromBits(uint_n(1)); }
  static constexpr dfloat infinity()      noexcept { return floatFromBits(uint_n( allOne - mantAllOne() )); }
  static constexpr dfloat quiet_NaN()     noexcept { return floatFromBits(uint_n( allOne - mantAllOne(-1) )); }
  static constexpr dfloat signaling_NaN() noexcept { return floatFromBits(uint_n( allOne - mantAllOne(-1) - mantAllOne(-2) )); }
};
}

#endif // DYNAMICFLOAT_H
