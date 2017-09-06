#ifndef DYNAMICFLOAT_H
#define DYNAMICFLOAT_H

#include <cstdint>
#include <iostream>
#include <limits>
#include <cstring>

template <typename TDest, typename TSource>
inline TDest bit_cast(const TSource& source)
{
  static_assert(sizeof(TDest) == sizeof(TSource),"Dest and source must have same size!");
  TDest dest;
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

template <unsigned TSignificand, unsigned TExponent>
class dynamicFloat
{
  //Friend the numeric_limits specialisation
  friend struct std::numeric_limits<dynamicFloat<TSignificand, TExponent>>;
  //Friend all specialisations of dynamic float for easier conversion
  template<unsigned, unsigned> friend class dynamicFloat;
  //Declare total number of bits for easy reference
  static constexpr int typeWidth = TSignificand + TExponent + 1;
  //Define the smallest unsigned int that can represent this
  using uint_n = uintN_t<typeWidth>;

public:

  dynamicFloat() = default;

  dynamicFloat(const float _inFloat) : dynamicFloat(bit_cast<dynamicFloat<23,8>>(_inFloat)){}
  dynamicFloat<TSignificand, TExponent>& operator= (const float _inFloat);

  dynamicFloat(const double _inDouble) : dynamicFloat(bit_cast<dynamicFloat<52,11>>(_inDouble)){}
  dynamicFloat<TSignificand, TExponent>& operator= (const double _inDouble);

  dynamicFloat(const dynamicFloat<TSignificand, TExponent> &_inFloat) : m_data(_inFloat.m_data.m_bits) {}
  dynamicFloat<TSignificand, TExponent>& operator= (const dynamicFloat<TSignificand, TExponent> &_inFloat);

  template<unsigned TInSigLen, unsigned TInExpLen>
  dynamicFloat(const dynamicFloat<TInSigLen,TInExpLen> _inFloat);

  operator float() const;

  operator double() const;

  bool operator== (const dynamicFloat<TSignificand, TExponent> _inFloat) const;

  template<unsigned TASigLen, unsigned TAExpLen, unsigned TBSigLen, unsigned TBExpLen>
  inline friend bool operator== (const dynamicFloat<TASigLen, TAExpLen> _a, const dynamicFloat<TBSigLen, TBExpLen> _b);

  template<unsigned TInSigLen, unsigned TInExpLen>
  inline friend bool operator== (const dynamicFloat<TInSigLen, TInExpLen> _a, const float _b);

  template<unsigned TInSigLen, unsigned TInExpLen>
  inline friend bool operator== (const float _a, const dynamicFloat<TInSigLen, TInExpLen> _b);

  template<unsigned TInSigLen, unsigned TInExpLen>
  inline friend bool operator== (const dynamicFloat<TInSigLen, TInExpLen> _a, const double _b);

  template<unsigned TInSigLen, unsigned TInExpLen>
  inline friend bool operator== (const double _a, const dynamicFloat<TInSigLen, TInExpLen> _b);

  template<unsigned TInSigLen, unsigned TInExpLen>
  auto operator* (const dynamicFloat<TInSigLen, TInExpLen> _rhs) const;

  template<unsigned TInSigLen, unsigned TInExpLen>
  auto operator/ (const dynamicFloat<TInSigLen, TInExpLen> _rhs) const;

  template<unsigned TInSigLen, unsigned TInExpLen>
  auto operator+ (const dynamicFloat<TInSigLen, TInExpLen> _rhs) const;

  template<unsigned TInSigLen, unsigned TInExpLen>
  auto operator- (const dynamicFloat<TInSigLen, TInExpLen> _rhs) const;

  dynamicFloat<TSignificand, TExponent> operator-() const;

  template<unsigned TInSigLen, unsigned TInExpLen>
  dynamicFloat<TSignificand, TExponent>& operator += (const dynamicFloat<TInSigLen, TInExpLen> _rhs);

  template<unsigned TInSigLen, unsigned TInExpLen>
  dynamicFloat<TSignificand, TExponent>& operator -= (const dynamicFloat<TInSigLen, TInExpLen> _rhs);

  //  template<unsigned TInSigLen, unsigned TInExpLen>
  //	dynamicFloat<TSignificand, TExponent>& operator *= (const dynamicFloat<TInSigLen, TInExpLen> _rhs);

  //  template<unsigned TInSigLen, unsigned TInExpLen>
  //	dynamicFloat<TSignificand, TExponent>& operator /= (const dynamicFloat<TInSigLen, TInExpLen> _rhs);

  //	dynamicFloat<TSignificand, TExponent>& operator += (float other);
  //	dynamicFloat<TSignificand, TExponent>& operator -= (float other);
  //	dynamicFloat<TSignificand, TExponent>& operator *= (float other);
  //	dynamicFloat<TSignificand, TExponent>& operator /= (float other);

  dynamicFloat<TSignificand, TExponent> operator--(const int) const;

  dynamicFloat<TSignificand, TExponent>& operator--() const;


private:
  // union to store this classes data
  union data
  {
    data() = default;
    data(uint_n _bits) : m_bits(_bits) {}
    uint_n m_bits;  // All bits
    struct
    {
      uint_n significand : TSignificand;	// mantissa
      uint_n exponent : TExponent;        // exponent
      uint_n sign : 1;                    // sign always 1 bit
    } IEEE;
  } m_data;

  static constexpr auto bias(const unsigned _exponent) noexcept;

  static constexpr auto minExponent() noexcept;

  static constexpr auto maxExponent(const unsigned _exponent) noexcept;

  static constexpr auto iShiftPow2(const unsigned _mantissa) noexcept;

  static constexpr unsigned msb(uint64_t v);

  static constexpr unsigned lsb(const uint64_t v);

  inline bool isNaN() const noexcept;

  inline bool isInfinity() const noexcept;

  inline bool isDenorm() const noexcept;

  template<int width>
  bool nonOverflowMult(const uint64_t _lhs, const uint64_t _rhs, uint64_t &io_result) const;

  template<int width>
  bool longDivide(uint64_t _lhs, const uint64_t _rhs, uint64_t &io_result) const;

  inline auto roundingShift(uint64_t io_num, const uint64_t _shift) const;

};


//----------------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------Numeric Limits Specialization------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------------

namespace std
{
template <unsigned TSignificand, unsigned TExponent>
struct numeric_limits<dynamicFloat<TSignificand, TExponent>>
{

private:
  using dfloat = dynamicFloat<TSignificand, TExponent>;
  static constexpr int typeWidthMinusOne = TSignificand + TExponent;
  static constexpr int typeWidth   = typeWidthMinusOne + 1;
  using uint_n = uintN_t<typeWidth>;
  static constexpr uint_n maskOnes = ~uint_n{0};
  static constexpr uint_n maskOnesNoSign = maskOnes & ~(1ul << typeWidthMinusOne); // 0-111.... , 0 followed by all ones
  static constexpr uint64_t mantAllOne(int64_t _offset) { return ((1ul << (TSignificand + _offset)) - 1 ); }
  static constexpr dfloat floatFromBits(uint_n _bits)
  {
    //Default constructor will avoid non-constexpr call to bit_cast
    dfloat bitFloat{};
    //Set the bits directly with friend class relationship
    bitFloat.m_data.m_bits = _bits;
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
  static constexpr int  digits       = TSignificand + 1;
  static constexpr int  digits10     = TSignificand * log_base10_radix;
  static constexpr int  max_digits10 = digits * log_base10_radix + 2;
  static constexpr bool is_signed		 = true;
  static constexpr bool is_integer	 = false;
  static constexpr bool is_exact		 = false;
  static constexpr bool traps        = false;
  static constexpr bool is_modulo		 = false;
  static constexpr bool is_bounded	 = true;

  // Floating point specific.

  static constexpr dfloat epsilon()     noexcept { return dfloat(1.0f / (1ul << TSignificand)); }
  static constexpr dfloat round_error() noexcept { return dfloat(0.5); }

  static constexpr int max_exponent   = 1 << (TExponent - 1);
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
  static constexpr dfloat signaling_NaN() noexcept { return floatFromBits(uint_n{maskOnesNoSign - mantAllOne(0) + mantAllOne(-2) + 1}); }
};
}

#include "./dynamicfloat.inl"

#endif // DYNAMICFLOAT_H
