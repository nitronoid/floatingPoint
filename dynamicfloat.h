#ifndef DYNAMICFLOAT_H
#define DYNAMICFLOAT_H

#include <cstdint>
#include <cmath>
#include <iostream>
#include <limits>

// Selects the correct size type based on passed bit num
template<int N>
using uintN_t =
  typename std::conditional_t< N == 8, uint8_t,
    typename std::conditional_t< N == 16, uint16_t,
      typename std::conditional_t< N == 32, uint32_t,
       uint64_t
    >
  >
>;

template <int MANTISSA, int EXPONENT>
class dynamicFloat
{
  static constexpr int DIGITS = MANTISSA + EXPONENT + 1;
  template<typename T>
  using uintT_t =
    typename std::conditional_t< std::is_same<T,uint8_t>::value, uint8_t,
      typename std::conditional_t< std::is_same<T,uint16_t>::value, uint16_t,
        typename std::conditional_t< std::is_same<T,float>::value || std::is_same<T,uint32_t>::value, uint32_t,
          uint64_t
      >
    >
  >;

  template<typename T>
  using df_t =
    typename std::conditional_t< std::is_same<T,uint8_t>::value, dynamicFloat<3,4>,
      typename std::conditional_t< std::is_same<T,uint16_t>::value, dynamicFloat<10,5>,
         typename std::conditional_t< std::is_same<T,float>::value || std::is_same<T,uint32_t>::value, dynamicFloat<23,8>,
          dynamicFloat<52,11>
      >
    >
  >;

public:

  dynamicFloat() = default;

  template<typename FPT>
  dynamicFloat(FPT _float32)
  {
    // define the IEEE precision type based on passed type
    // get the mantissa length
    constexpr int ms = std::numeric_limits<df_t<FPT>>::digits - 1;
    // get the exponent length
    constexpr int es = sizeof(FPT)*8 - ms - 1;
    // define an alias for easy reference to the type
    using IEEE_t = IEEET<FPT, uintT_t<FPT>, ms, es>;
    // construct the type from passed in variable
    IEEE_t floatN(_float32);

    // sign bit should remain the same
    IEEE.sign = floatN.IEEE.sign;

    // calculate the exponent value without bias, single has 127 bias
    int64_t unbiasedExponent = floatN.IEEE.exponent-bias(IEEE_t::expLen());

    // find the difference in size between the mantissas
    int64_t padSize = padLen<IEEE_t>();

    if(unbiasedExponent < minExponent())
    {
      IEEE.mantissa = IEEE.exponent = 0;
    }

    // if NaN or Infinity, stay NaN or Infinity
    else if (floatN.IEEE.exponent == maxExp(IEEE_t::expLen()))
    {
      // NaN or INF
      IEEE.mantissa = (floatN.IEEE.mantissa != 0) ? 1 : 0;
      IEEE.exponent = maxExp(EXPONENT);
    }

    else if (unbiasedExponent<1-bias(EXPONENT))
    {
      // this maps to a subnormal number
      IEEE.exponent = 0;
      uint64_t shiftVal =  (unbiasedExponent + bias(EXPONENT));
      IEEE.mantissa = (base(MANTISSA) >> shiftVal) + ( floatN.IEEE.mantissa >> (padSize + shiftVal));
    }

    // if the exponent is zero, map to zero
    else if (!floatN.IEEE.exponent)
    {
      if(!floatN.IEEE.mantissa)
        IEEE.mantissa = IEEE.exponent = 0;
      else
      {
        int64_t e = 0;
        uint64_t m = floatN.IEEE.mantissa<<1;

        for(;!(m & base(IEEE_t::mantLen())); ++e, m <<= 1);

        IEEE.exponent = bias(EXPONENT) - bias(IEEE_t::expLen()) - e;
        IEEE.mantissa = (padSize < 0) ? ((m & (base(IEEE_t::mantLen())-1)) << -padSize ) : ((m & (base(IEEE_t::mantLen())-1)) >> padSize);
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
      uint64_t foo = floatN.IEEE.mantissa;
      IEEE.mantissa = (padSize < 0) ? (foo << -padSize) : (foo >> padSize);
    }
  }

  operator float() const
  {
    dynamicFloat<23,8> f(m_bits);
    uint32_t foo = f.getBits();
    return *((float*)&foo);
  }

  operator double() const
  {
    dynamicFloat<52,11> f(m_bits);
    uint64_t foo = f.getBits();
    return *((double*)&foo);
  }

  inline uintN_t<DIGITS> getBits() const noexcept
  {
    return m_bits;
  }

private:
  // union to store this classes data
  union
  {
    uintN_t<DIGITS> m_bits;                // All bits
    struct
    {
      uintN_t<DIGITS> mantissa : MANTISSA;	// mantissa
      uintN_t<DIGITS> exponent : EXPONENT;	// exponent
      uintN_t<DIGITS> sign : 1;            // sign always 1 bit
    } IEEE;
  };

  template<typename data_t, typename raw_t,int mant, int exp>
  union IEEET
  {
    IEEET() = default;
    IEEET(data_t _data) : data(_data) {}
    IEEET(raw_t _m, raw_t _e, raw_t _s) {IEEE.mantissa = _m; IEEE.exponent = _e; IEEE.sign = _s;}
    data_t data;
    struct
    {
      raw_t mantissa : mant;
      raw_t exponent  : exp;
      raw_t sign : 1;
    } IEEE;
    static constexpr int mantLen() noexcept { return mant; }
    static constexpr int expLen() noexcept { return exp; }
  };

  static constexpr int64_t bias(uint64_t _exponent) noexcept { return (uint64_t(1) << (_exponent-1)) - 1; }

  static constexpr int64_t minExponent() noexcept { return -(MANTISSA - 1 + bias(EXPONENT)); }

  static constexpr uint64_t maxExp(uint64_t _exponent) noexcept { return (uint64_t(1) << _exponent) - 1; }

  static constexpr uint64_t base(uint64_t _mantissa) noexcept { return uint64_t(1) << _mantissa; }

  template<typename T>
  static constexpr int64_t padLen() noexcept { return T::mantLen() - MANTISSA; }
};

namespace std
{
template <int MANTISSA, int EXPONENT>
class numeric_limits<dynamicFloat<MANTISSA, EXPONENT>>
{
  using dfloat = dynamicFloat<MANTISSA, EXPONENT>;
  static constexpr int calcDigits10() noexcept { return MANTISSA * 0.30102999566; } // *log10(radix)
  static constexpr int N = MANTISSA + EXPONENT + 1;
public:

  // General -- meaningful for all specializations.

  static const bool is_specialized = true;
  static constexpr dfloat min () { return dfloat(uintN_t<N>(uint64_t(1)<<MANTISSA)); }
  static constexpr dfloat max () {return dfloat(uintN_t<N>((uint64_t(1)<<(N-1))-1-(uint64_t(1)<<MANTISSA)));}
  static const int radix = 2;
  static const int digits = MANTISSA + 1;   // conservative assumption
  static const int digits10 = calcDigits10();  // conservative assumption
  static const bool is_signed		= true;
  static const bool is_integer	= false;
  static const bool is_exact		= false;
  static const bool traps       = false;
  static const bool is_modulo		= false;
  static const bool is_bounded	= true;

  // Floating point specific.

  static constexpr dfloat epsilon () noexcept { return dfloat(pow(2,-MANTISSA)); } //Inefficient!
  static constexpr dfloat round_error () noexcept {return dfloat(0.5);}
  //   static const int min_exponent10 = HalfFloat::MIN_EXPONENT10;
  //   static const int max_exponent10 = HalfFloat::MAX_EXPONENT10;
  static constexpr int min_exponent   = 1 - ((1 << (EXPONENT-1))-2);
  static const int max_exponent   = 1<<(EXPONENT-1);

  static const bool has_infinity			= true;
  static const bool has_quiet_NaN			= true;
  static const bool has_signaling_NaN		= true;
  static const bool is_iec559				= false;
  static const bool has_denorm			= denorm_present;
  static const bool tinyness_before		= false;
  static const float_round_style round_style = round_to_nearest;

  static constexpr dfloat denorm_min () {return dfloat(uintN_t<N>(1));}
  static dfloat infinity () {return dfloat(uintN_t<N>((uint64_t(1)<<(N-1))-1 - ((uint64_t(1)<<(MANTISSA))-1)));}
  //   static HalfFloat quiet_NaN ()
  //   {return HalfFloat(1,HalfFloat::MAX_EXPONENT_VALUE,0);}
  //   static HalfFloat signaling_NaN ()
  //   {return HalfFloat(1,HalfFloat::MAX_EXPONENT_VALUE,0);}
};
}

#endif // DYNAMICFLOAT_H
