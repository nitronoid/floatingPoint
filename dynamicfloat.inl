
#ifndef DYNAMICFLOAT_INL_INCLUDED
#define DYNAMICFLOAT_INL_INCLUDED


//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
dynamicFloat<TSignificand, TExponent>& dynamicFloat<TSignificand, TExponent>::operator= (const float _inFloat)
{
  m_bits = dynamicFloat<TSignificand, TExponent>(_inFloat).m_bits;
  return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
dynamicFloat<TSignificand, TExponent>& dynamicFloat<TSignificand, TExponent>::operator= (const double _inDouble)
{
  m_bits = dynamicFloat<TSignificand, TExponent>(_inDouble).m_bits;
  return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
dynamicFloat<TSignificand, TExponent>& dynamicFloat<TSignificand, TExponent>::operator= (const dynamicFloat<TSignificand, TExponent> &_inFloat)
{
  m_bits = _inFloat.m_bits;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
template<unsigned TInSigLen, unsigned TInExpLen>
dynamicFloat<TSignificand, TExponent>::dynamicFloat(const dynamicFloat<TInSigLen,TInExpLen> _inFloat)
{
  // sign bit should remain the same
  IEEE.sign = _inFloat.IEEE.sign;

  // calculate the exponent value without bias, single has 127 bias
  auto unbiasedExponent = _inFloat.IEEE.exponent - bias(TInExpLen);

  // find the difference in size between the mantissas
  constexpr auto padSize = static_cast<int>(TInSigLen - TSignificand);

  // too small to represent
  if(unbiasedExponent < minExponent())
  {
    IEEE.significand = IEEE.exponent = 0;
  }

  // if NaN or Infinity, stay NaN or Infinity
  else if (_inFloat.IEEE.exponent == maxExponent(TInExpLen))
  {
    IEEE.significand = _inFloat.IEEE.significand;
    IEEE.exponent = maxExponent(TExponent);
  }

  // this maps to a subnormal number
  else if (unbiasedExponent < (1 - bias(TExponent)))
  {
    IEEE.exponent = 0;
    auto shiftVal = 1 - (unbiasedExponent + bias(TExponent));
    IEEE.significand = (iShiftPow2(TSignificand) >> shiftVal) + ( _inFloat.IEEE.significand >> (padSize + shiftVal));
  }

  // if the exponent is zero, map to zero
  else if (!_inFloat.IEEE.exponent)
  {
    if(!_inFloat.IEEE.significand)
      IEEE.significand = IEEE.exponent = 0;
    else
    {
      auto base = iShiftPow2(TInSigLen);
      auto e = TInSigLen - msb(_inFloat.IEEE.significand);
      auto m = _inFloat.IEEE.significand << e;
      IEEE.exponent = bias(TExponent) - bias(TInExpLen) - e + 1;
      IEEE.significand = (padSize < 0) ?
            ((m & (base - 1)) << -padSize)
          : ((m & (base - 1)) >> padSize);
    }
  }

  // too large to be represented, map this value to infinity
  else if (unbiasedExponent > bias(TExponent))
  {
    IEEE.significand = 0;
    IEEE.exponent = maxExponent(TExponent);
  }

  // normal numbers
  else
  {
    // add the new bias
    IEEE.exponent = unbiasedExponent + bias(TExponent);
    auto signicandCopy = static_cast<uint64_t>(_inFloat.IEEE.significand);
    IEEE.significand = (padSize < 0) ?
          (signicandCopy << -padSize)
        : (signicandCopy >> padSize);
  }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
dynamicFloat<TSignificand, TExponent>::operator float() const
{
  return bit_cast<float>(dynamicFloat<23,8>(*this));
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
dynamicFloat<TSignificand, TExponent>::operator double() const
{
  return bit_cast<double>(dynamicFloat<52,11>(*this));
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
bool dynamicFloat<TSignificand, TExponent>::operator== (const dynamicFloat<TSignificand, TExponent> _inFloat) const
{
  // +0 and -0 are considered to be equal
  if (!(m_bits << 1u) && !(_inFloat.m_bits << 1u)) return true;
  return m_bits == _inFloat.m_bits && !isNaN();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<unsigned TASigLen, unsigned TAExpLen, unsigned TBSigLen, unsigned TBExpLen>
bool operator== (const dynamicFloat<TASigLen, TAExpLen> _a, const dynamicFloat<TBSigLen, TBExpLen> _b)
{
  dynamicFloat<TASigLen, TAExpLen> converted = _b;
  return _a == converted;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<unsigned TInSigLen, unsigned TInExpLen>
bool operator== (const dynamicFloat<TInSigLen, TInExpLen> _a, const float _b)
{
  return _a == dynamicFloat<TInSigLen, TInExpLen>(_b);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<unsigned TInSigLen, unsigned TInExpLen>
bool operator== (const float _a, const dynamicFloat<TInSigLen, TInExpLen> _b)
{
  return _b == _a;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<unsigned TInSigLen, unsigned TInExpLen>
bool operator== (const dynamicFloat<TInSigLen, TInExpLen> _a, const double _b)
{
  return _a == dynamicFloat<TInSigLen, TInExpLen>(_b);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<unsigned TInSigLen, unsigned TInExpLen>
bool operator== (const double _a, const dynamicFloat<TInSigLen, TInExpLen> _b)
{
  return _b == _a;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
template<unsigned TInSigLen, unsigned TInExpLen>
auto dynamicFloat<TSignificand, TExponent>::operator+ (const dynamicFloat<TInSigLen, TInExpLen> _lhs)
{
  // define the return type as the most precision passed
  using TLargest = typename std::conditional_t< ((TSignificand + TExponent) > (TInSigLen + TInExpLen)), dynamicFloat<TSignificand, TExponent>, dynamicFloat<TInSigLen, TInExpLen>>;
  constexpr unsigned siglen = std::max(TSignificand,TInSigLen);
  TLargest rhsConv = *this;
  TLargest lhsConv { _lhs};

  //if either is zero return the other
  if(!rhsConv.m_bits) return lhsConv;
  if(!lhsConv.m_bits) return rhsConv;

  // if either passed is NaN the result is NaN
  if(rhsConv.isNaN() || lhsConv.isNaN())
  {
    return std::numeric_limits<TLargest>::signaling_NaN();
  }

  // if either passed is infinity the result is infinity
  else if(rhsConv.isInfinity() || lhsConv.isInfinity())
  {
    return std::numeric_limits<TLargest>::infinity();
  }

  TLargest result;
  // calculate the difference between exponents
  int expDiff = rhsConv.IEEE.exponent - lhsConv.IEEE.exponent;

  // get the significands and add the hidden bit,
  // IEEE uses 1.xxx... where the 1 isn't stored, so we add it
  auto rhsSig = rhsConv.IEEE.significand | (uint64_t{1} << siglen);
  auto lhsSig = lhsConv.IEEE.significand | (uint64_t{1} << siglen);

  // get both operands in the same format with same exponent
  if(expDiff > 0)
  {
    lhsSig >>= expDiff;
    // choose this format's exponent
    result.IEEE.exponent = rhsConv.IEEE.exponent;
  }
  else
  {
    rhsSig >>= -expDiff;
    // choose this format's exponent
    result.IEEE.exponent = lhsConv.IEEE.exponent;
  }

  // use two's complement for negative
  if(rhsConv.IEEE.sign) rhsSig = -rhsSig;
  if(lhsConv.IEEE.sign) lhsSig = -lhsSig;

  // sum the significands for our result
  auto resultSig = rhsSig + lhsSig;

  // path dependent on sign of result
  if(resultSig >> (sizeof(TLargest) * 8 - 1))
  {
    // negative sign
    result.IEEE.sign = 1;
    // convert result significand to two's complement
    result.IEEE.significand = -resultSig;
  }
  else
  {
    // amount of overflow bits we need to normalize
    auto overflow = msb(resultSig) - siglen;
    // check sign of new significand
    if(overflow >> (sizeof(overflow) * 8 - 1))
      // shift with rounding
      resultSig = (resultSig >> overflow) + msb((resultSig & (uint64_t{1}<<overflow)) + 1);
    else
      // don't need to round, reverse shift for negative
      resultSig <<= -overflow;

    // positive sign
    result.IEEE.sign = 0;
    // assign our calculated significand
    result.IEEE.significand = resultSig;
    // add the overflow to our exponent to normalize the float
    result.IEEE.exponent += overflow;
  }
  return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
bool dynamicFloat<TSignificand, TExponent>::isNaN() const
{
  return IEEE.significand != 0 && IEEE.exponent == maxExponent(TExponent);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
bool dynamicFloat<TSignificand, TExponent>::isInfinity() const
{
  return IEEE.significand == 0 && IEEE.exponent == maxExponent(TExponent);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
bool dynamicFloat<TSignificand, TExponent>::isDenorm() const
{
  return IEEE.exponent == 0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
constexpr auto dynamicFloat<TSignificand, TExponent>::bias(unsigned _exponent) noexcept
{
  return static_cast<int>(1ul << (_exponent - 1)) - 1;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
constexpr auto dynamicFloat<TSignificand, TExponent>::minExponent() noexcept
{
  return -static_cast<int>(TSignificand - 1 + bias(TExponent));
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
constexpr auto dynamicFloat<TSignificand, TExponent>::maxExponent(unsigned _exponent) noexcept
{
  return (1ul << _exponent) - 1;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
constexpr auto dynamicFloat<TSignificand, TExponent>::iShiftPow2(unsigned _mantissa) noexcept
{
  return 1ul << _mantissa;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
constexpr auto dynamicFloat<TSignificand, TExponent>::ilog2(unsigned n, unsigned p) noexcept
{
  return (n <= 1) ? p : ilog2(n / 2, p + 1);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
constexpr uint32_t dynamicFloat<TSignificand, TExponent>::msb(uint32_t v)
{
  constexpr int deBruijnTable[32] = {
    0,   9,   1,  10,  13,  21,   2,  29,
    11,  14,  16,  18,  22,  25,   3,  30,
    8,  12,  20,  28,  15,  17,  24,   7,
    19,  27,  23,   6,  26,   5,   4,  31
  };

  // is constant time
  for(int i = 1; i < 32; i <<= 1) v |= v >> i;

  return deBruijnTable[static_cast<uint32_t>( v * 0x07C4ACDDU ) >> 27];
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#endif
