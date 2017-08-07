#ifndef DYNAMICFLOAT_H
#define DYNAMICFLOAT_H

#include <cstdint>
#include <cmath>
#include <iostream>
#include <limits>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wbitfield-constant-conversion"

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

template<typename T>
using uintT_t =
    typename std::conditional_t< std::is_same<T,uint8_t>::value, uint8_t,
        typename std::conditional_t< std::is_same<T,uint16_t>::value, uint16_t,
            typename std::conditional_t< std::is_same<T,float>::value || std::is_same<T,uint32_t>::value, uint32_t,
                uint64_t
            >
        >
    >;

template <int MANTISSA, int EXPONENT, int N = MANTISSA + EXPONENT + 1>
class dynamicFloat
{

public:

  dynamicFloat() = default;

  template<typename FPT>
  dynamicFloat(FPT _float32)
  {
    // create a single precision float
    //IEEEn<FPT> floatN(_float32);
    constexpr uint64_t ms = std::numeric_limits<float>::digits - 1;
    constexpr uint64_t es = sizeof(float) - ms - 1;
    using IEEE_t = IEEET<FPT, uintT_t<FPT>, ms, es>;
    IEEE_t floatN(_float32);

    // sign bit should remain the same
    IEEE.sign = floatN.IEEE.sign;

    // calculate the exponent value without bias, single has 127 bias
    int64_t unbiasedExponent = floatN.IEEE.exponent-bias(expLen<IEEE_t>());

    // find the difference in size between the mantissas
    int64_t padSize = padLen<IEEE_t>();

    // if the exponent is zero, map to zero
    if (!floatN.IEEE.exponent)
    {
      if(!floatN.IEEE.mantissa)
        IEEE.mantissa = IEEE.exponent = 0;
      else
      {
        int64_t e = 0;
        uint64_t m = floatN.IEEE.mantissa<<1;


        for(;!(m & base(mantLen<IEEE_t>())); ++e, m <<= 1);

        IEEE.exponent = bias(EXPONENT) - bias(expLen<IEEE_t>()) - e;
        IEEE.mantissa = (padSize < 0) ? ((m & (base(mantLen<IEEE_t>())-1)) << -padSize ) : ((m & (base(mantLen<IEEE_t>())-1)) >> padSize);
      }
    }

    else if(unbiasedExponent < minExponent())
    {
      IEEE.mantissa = IEEE.exponent = 0;
    }

    // if NaN or Infinity, stay NaN or Infinity
    else if (floatN.IEEE.exponent == maxExp(expLen<IEEE_t>()))
    {
      // NaN or INF
      IEEE.mantissa = (floatN.IEEE.mantissa != 0) ? 1 : 0;
      IEEE.exponent = maxExp(EXPONENT);
    }

    else if (unbiasedExponent<1-bias(EXPONENT))
    {
      // this maps to a subnormal number
      IEEE.exponent = 0;
      uint64_t shiftVal =  (1-bias(EXPONENT) - unbiasedExponent);  // 2^-exp_val
      IEEE.mantissa = (base(MANTISSA) >> shiftVal) + ( floatN.IEEE.mantissa >> (padSize + shiftVal));
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

  inline uintN_t<N> getBits() const noexcept
  {
    return m_bits;
  }

private:
  // union to store this classes data
  union
  {
    uintN_t<N> m_bits;                // All bits
    struct
    {
      uintN_t<N> mantissa : MANTISSA;	// mantissa
      uintN_t<N> exponent : EXPONENT;	// exponent
      uintN_t<N> sign : 1;            // sign always 1 bit
    } IEEE;
  };

  template<typename f, typename d,int m, int e>
  union IEEET
  {
    IEEET() = default;
    IEEET(f _data) : data(_data) {}
    IEEET(d _m, d _e, d _s) {IEEE.mantissa = _m; IEEE.exponent = _e; IEEE.sign = _s;}
    f data;
    struct
    {
      d mantissa : m;
      d exponent  : e;
      d sign : 1;
    } IEEE;
  };

  union IEEEMini
  {
    IEEEMini() = default;
    IEEEMini(uint8_t _data) : data(_data) {}
    IEEEMini(uint8_t _m, uint8_t _e, uint8_t _s) {IEEE.mantissa = _m; IEEE.exponent = _e; IEEE.sign = _s;}
    uint8_t data;
    struct
    {
      uint8_t mantissa : 3;
      uint8_t exponent  : 4;
      uint8_t sign : 1;
    } IEEE;
  };

  union IEEEHalf
  {
    IEEEHalf() = default;
    IEEEHalf(uint16_t _data) : data(_data) {}
    IEEEHalf(uint16_t _m, uint16_t _e, uint16_t _s) {IEEE.mantissa = _m; IEEE.exponent = _e; IEEE.sign = _s;}
    uint16_t data;
    struct
    {
      uint16_t mantissa : 10;
      uint16_t exponent  : 5;
      uint16_t sign : 1;
    } IEEE;
  };

  // union to access float bits
  union IEEESingle
  {
    IEEESingle() = default;
    IEEESingle(float _float) : float_32(_float) {}
    IEEESingle(uint32_t _data) : data(_data) {}
    IEEESingle(uint32_t _m, uint32_t _e, uint32_t _s) {IEEE.mantissa = _m; IEEE.exponent = _e; IEEE.sign = _s;}
    float float_32;
    uint32_t data;
    struct
    {
      uint32_t mantissa : 23;
      uint32_t exponent  : 8;
      uint32_t sign : 1;
    } IEEE;
  };

  // union to access double bits
  union IEEEDouble
  {
    IEEEDouble() = default;
    IEEEDouble(double _float) : double_64(_float) {}
    IEEEDouble(uint64_t _data) : data(_data) {}
    IEEEDouble(uint64_t _m, uint64_t _e, uint64_t _s) {IEEE.mantissa = _m; IEEE.exponent = _e; IEEE.sign = _s;}
    double double_64;
    uint64_t data;
    struct
    {
      uint64_t mantissa : 52;
      uint64_t exponent  : 11;
      uint64_t sign : 1;
    } IEEE;
  };

  template<typename T>
  using IEEEn =
    typename std::conditional_t< std::is_same<T,float>::value, IEEESingle,
        typename std::conditional_t< std::is_same<T,double>::value, IEEEDouble,
            typename std::conditional_t< std::is_same<T,uint32_t>::value, IEEESingle,
                typename std::conditional_t< std::is_same<T,uint64_t>::value, IEEEDouble,
                    typename std::conditional_t< std::is_same<T,uint16_t>::value, IEEEHalf,
                        IEEEMini
                    >
                >
            >
        >
    >;

  static constexpr int64_t ipow(int64_t base, int64_t exp, int64_t result = 1) noexcept
  {
    return exp < 1 ? result : ipow(base*base, exp/2, (exp % 2) ? result*base : result);
  }

  static constexpr int64_t bias(uint64_t _exponent) noexcept { return ipow(2, _exponent - 1) - 1; }

  static constexpr uint64_t bitLen(uint64_t _bitField) noexcept
  {
    unsigned int len = 0;
    for(;_bitField; _bitField &= (_bitField - 1), ++len);
    return len;
  }

  template<typename T>
  uint64_t expLen() noexcept
  {
    T _bitField;
    _bitField.IEEE.exponent = ~0;
    return bitLen(_bitField.IEEE.exponent);
  }

  template<typename T>
  uint64_t mantLen() noexcept
  {
    T _bitField;
    _bitField.IEEE.mantissa = ~0;
    return bitLen(_bitField.IEEE.mantissa);
  }

  static constexpr int64_t minExponent() noexcept { return -(MANTISSA - 1 + bias(EXPONENT)); }

  static constexpr uint64_t maxExp(uint64_t _exponent) noexcept { return ipow(2, _exponent + 1) - 1; }

  static constexpr uint64_t base(uint64_t _mantissa) noexcept { return ipow(2, _mantissa); }

  template<typename T>
  int64_t padLen() const noexcept
  {
    T _bitField;
    _bitField.IEEE.mantissa = ~0;
    return bitLen(_bitField.IEEE.mantissa) - MANTISSA;
  }
};

#pragma clang diagnostic pop

#endif // DYNAMICFLOAT_H
