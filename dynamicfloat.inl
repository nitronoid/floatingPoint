
#ifndef DYNAMICFLOAT_INL_INCLUDED
#define DYNAMICFLOAT_INL_INCLUDED

//==================================================================================================================
//---------------------------------------------------Constructors---------------------------------------------------
//==================================================================================================================

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
dynamicFloat<TSignificand, TExponent>& dynamicFloat<TSignificand, TExponent>::operator= (const float _inFloat)
{
  m_data.m_bits = dynamicFloat<TSignificand, TExponent>(_inFloat).m_data.m_bits;
  return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
dynamicFloat<TSignificand, TExponent>& dynamicFloat<TSignificand, TExponent>::operator= (const double _inDouble)
{
  m_data.m_bits = dynamicFloat<TSignificand, TExponent>(_inDouble).m_data.m_bits;
  return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
dynamicFloat<TSignificand, TExponent>& dynamicFloat<TSignificand, TExponent>::operator= (const dynamicFloat<TSignificand, TExponent> &_inFloat)
{
  m_data.m_bits = _inFloat.m_data.m_bits;
  return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
template<unsigned TInSigLen, unsigned TInExpLen>
dynamicFloat<TSignificand, TExponent>::dynamicFloat(const dynamicFloat<TInSigLen,TInExpLen> _inFloat)
{
  // sign bit should remain the same
  m_data.IEEE.sign = _inFloat.m_data.IEEE.sign;

  // calculate the exponent value without bias, single has 127 bias
  auto unbiasedExponent = _inFloat.m_data.IEEE.exponent - bias(TInExpLen);

  // find the difference in size between the mantissas
  constexpr auto padSize = static_cast<int>(TInSigLen - TSignificand);

  // too small to represent
  if(unbiasedExponent < minExponent())
  {
    m_data.IEEE.significand = m_data.IEEE.exponent = 0;
  }

  // if NaN or Infinity, stay NaN or Infinity
  else if (_inFloat.m_data.IEEE.exponent == maxExponent(TInExpLen))
  {
    m_data.IEEE.significand = _inFloat.m_data.IEEE.significand;
    m_data.IEEE.exponent = maxExponent(TExponent);
  }

  // when the exponents are the same we need to renormalize the significand
  else if (TExponent == TInExpLen)
  {
    m_data.IEEE.exponent = _inFloat.m_data.IEEE.exponent;
    m_data.IEEE.significand = normalize<TSignificand>(_inFloat.m_data.IEEE.significand);
  }

  // this maps to a subnormal number
  else if (unbiasedExponent < (1 - bias(TExponent)))
  {
    m_data.IEEE.exponent = 0;
    auto shiftVal = 1 - (unbiasedExponent + bias(TExponent));
    m_data.IEEE.significand = (iShiftPow2(TSignificand) >> shiftVal) + ( _inFloat.m_data.IEEE.significand >> (padSize + shiftVal));
  }

  // if the exponent is zero, map to zero
  else if (!_inFloat.m_data.IEEE.exponent)
  {
    if(!_inFloat.m_data.IEEE.significand)
      m_data.IEEE.significand = m_data.IEEE.exponent = 0;
    else
    {
      auto base = iShiftPow2(TInSigLen);
      auto e = TInSigLen - (msb(_inFloat.m_data.IEEE.significand)-1);
      auto m = static_cast<unsigned>(_inFloat.m_data.IEEE.significand << e);
      m_data.IEEE.exponent = bias(TExponent) - bias(TInExpLen) - e + 1;
      m_data.IEEE.significand = (padSize < 0) ?
            ((m & (base - 1)) << -padSize)
          : ((m & (base - 1)) >> padSize);
    }
  }

  // too large to be represented, map this value to infinity
  else if (unbiasedExponent > bias(TExponent))
  {
    m_data.IEEE.significand = 0;
    m_data.IEEE.exponent = maxExponent(TExponent);
  }

  // normal numbers
  else
  {
    // add the new bias
    m_data.IEEE.exponent = static_cast<uint_n>(unbiasedExponent + bias(TExponent));
    auto signicandCopy = static_cast<uint64_t>(_inFloat.m_data.IEEE.significand);
    m_data.IEEE.significand = (padSize < 0) ?
          (signicandCopy << -padSize)
        : (signicandCopy >> padSize);
  }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//==================================================================================================================
//----------------------------------------------------Operators-----------------------------------------------------
//==================================================================================================================

//---------------------------------------------------Conversion-----------------------------------------------------

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
dynamicFloat<TSignificand, TExponent>::operator float() const
{
  return detail::bit_cast<float>(dynamicFloat<23,8>(*this));
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
dynamicFloat<TSignificand, TExponent>::operator double() const
{
  return detail::bit_cast<double>(dynamicFloat<52,11>(*this));
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------Equality------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
bool dynamicFloat<TSignificand, TExponent>::operator== (const dynamicFloat<TSignificand, TExponent> _inFloat) const
{
  // +0 and -0 are considered to be equal
  auto  thisNoSign = m_data.m_bits << 1u;
  auto inNoSign = _inFloat.m_data.m_bits << 1u;
  if (!thisNoSign && !inNoSign) return true;
  return m_data.m_bits == _inFloat.m_data.m_bits && !isNaN();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<unsigned TASigLen, unsigned TAExpLen, unsigned TBSigLen, unsigned TBExpLen>
bool operator== (const dynamicFloat<TASigLen, TAExpLen> _a, const dynamicFloat<TBSigLen, TBExpLen> _b)
{
  dynamicFloat<TASigLen, TAExpLen> converted = _b;
  return _a == converted;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------Arithmetic----------------------------------------------------

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
template<unsigned TInSigLen, unsigned TInExpLen>
auto dynamicFloat<TSignificand, TExponent>::operator* (const dynamicFloat<TInSigLen, TInExpLen> _rhs) const
{
  // define the return type as the most precision passed
  using TLargest = typename std::conditional_t<((TSignificand + TExponent) > (TInSigLen + TInExpLen)), dynamicFloat<TSignificand, TExponent>, dynamicFloat<TInSigLen, TInExpLen>>;
  // get the length of the significand and exponent for our resulting type
  constexpr unsigned siglen = std::max(TSignificand,TInSigLen);
  constexpr unsigned explen = std::max(TExponent,TInExpLen);
  // convert both operands to the result type
  TLargest lhsConv = *this;
  TLargest rhsConv = _rhs;
  TLargest result;

  // if either is zero return zero
  if(!lhsConv.m_data.m_bits || !rhsConv.m_data.m_bits) return TLargest{0.f};
  // if either passed is NaN or infinity return the same
  if(lhsConv.isNaN() || lhsConv.isInfinity()) return lhsConv;
  if(rhsConv.isNaN() || rhsConv.isInfinity()) return rhsConv;

  // Declare a variable of value one, with the same width as the result significand
  // This saves time during arithmetic operations, by not using the largest width for smaller types
  constexpr decltype(result.m_data.IEEE.significand) one = 1;

  // get the significands and add the hidden bit,
  // m_data.IEEE uses 1.xxx... where the 1 isn't stored, so we add it
  auto lhsSig = lhsConv.m_data.IEEE.significand | (one << siglen);
  auto rhsSig = rhsConv.m_data.IEEE.significand | (one << siglen);

  // exponent is the sum of the operands exponents
  result.m_data.IEEE.exponent = (lhsConv.m_data.IEEE.exponent - bias(explen)) + (rhsConv.m_data.IEEE.exponent - bias(explen)) + bias(explen);

  // make the same width as current signifcands
  // multiply for result
  decltype(lhsSig) resultSig;
  // returns whether the exponent needs to compensate for normalization
  bool up = nonOverflowMult<(siglen + 1) * 2>(lhsSig, rhsSig, resultSig);

  // XOR sign
  result.m_data.IEEE.sign = (lhsConv.m_data.IEEE.sign != rhsConv.m_data.IEEE.sign);
  // normalize the result
  result.m_data.IEEE.significand = normalize<siglen>(resultSig);
  // avoid a branch here and just add the bool
  result.m_data.IEEE.exponent += static_cast<unsigned>(up);
  return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
template<unsigned TInSigLen, unsigned TInExpLen>
auto dynamicFloat<TSignificand, TExponent>::operator/ (const dynamicFloat<TInSigLen, TInExpLen> _rhs) const
{
  // define the return type as the most precision passed
  using TLargest = typename std::conditional_t<((TSignificand + TExponent) > (TInSigLen + TInExpLen)), dynamicFloat<TSignificand, TExponent>, dynamicFloat<TInSigLen, TInExpLen>>;
  // get the length of the significand and exponent for our resulting type
  constexpr unsigned siglen = std::max(TSignificand,TInSigLen);
  constexpr unsigned explen = std::max(TExponent,TInExpLen);
  // convert both operands to the result type
  TLargest lhsConv = *this;
  TLargest rhsConv = _rhs;
  TLargest result;

  // if right is zero return NaN
  if(!rhsConv.m_data.m_bits) return std::numeric_limits<TLargest>::signaling_NaN();
  // if left is zero return zero
  if(!lhsConv.m_data.m_bits) return TLargest{0.0};
  // if either passed is NaN or infinity return the same
  if(lhsConv.isNaN() || lhsConv.isInfinity()) return lhsConv;
  if(rhsConv.isNaN() || rhsConv.isInfinity()) return rhsConv;

  // Declare a variable of value one, with the same width as the result significand
  // This saves time during arithmetic operations, by not using the largest width for smaller types
  constexpr decltype(result.m_data.IEEE.significand) one = 1;

  // get the significands and add the hidden bit,
  // m_data.IEEE uses 1.xxx... where the 1 isn't stored, so we add it
  auto lhsSig = lhsConv.m_data.IEEE.significand | (one << siglen);
  auto rhsSig = rhsConv.m_data.IEEE.significand | (one << siglen);

  // exponent is the difference of the operands exponents
  result.m_data.IEEE.exponent = (lhsConv.m_data.IEEE.exponent - bias(explen)) - (rhsConv.m_data.IEEE.exponent - bias(explen)) + bias(explen);
  // divide the significands for our result
  decltype(lhsSig) resultSig;
  // returns whether the exponent needs to compensate for normalization
  bool down = longDivide<siglen + 1>(lhsSig , rhsSig, resultSig);

  // XOR sign
  result.m_data.IEEE.sign = (lhsConv.m_data.IEEE.sign != rhsConv.m_data.IEEE.sign);
  // normalize the result
  result.m_data.IEEE.significand = normalize<siglen>(resultSig);
  // avoid a branch here and just add the bool
  result.m_data.IEEE.exponent -= static_cast<unsigned>(down);
  return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
template<unsigned TInSigLen, unsigned TInExpLen>
auto dynamicFloat<TSignificand, TExponent>::operator+ (const dynamicFloat<TInSigLen, TInExpLen> _rhs) const
{
  // define the return type as the most precision passed
  using TLargest = typename std::conditional_t<((TSignificand + TExponent) > (TInSigLen + TInExpLen)), dynamicFloat<TSignificand, TExponent>, dynamicFloat<TInSigLen, TInExpLen>>;
  // get the length of the significand and exponent for our resulting type
  constexpr unsigned siglen = std::max(TSignificand,TInSigLen);
  constexpr unsigned explen = std::max(TExponent,TInExpLen);
  // convert both operands to the result type
  TLargest lhsConv = *this;
  TLargest rhsConv = _rhs;
  TLargest result;

  // if left side is infinity, NaN, or right is zero, return the left
  if(lhsConv.isNaN() || lhsConv.isInfinity() || !rhsConv.m_data.m_bits) return lhsConv;
  // if right side is infinity, NaN, or right is zero, return the right
  if(rhsConv.isNaN() || rhsConv.isInfinity() || !lhsConv.m_data.m_bits) return rhsConv;

  // Declare a variable of value one, with the same width as the result significand
  // This saves time during arithmetic operations, by not using the largest width for smaller types
  constexpr decltype(result.m_data.IEEE.significand) one = 1;
  constexpr int signMap[2] = {1,-1};

  // get the significands and add the hidden bit,
  // m_data.IEEE uses 1.xxx... where the 1 isn't stored, so we add it
  // we also bias the result to avoid precision loss in the next step
  auto lhsSig = (lhsConv.m_data.IEEE.significand | (one << siglen)) << (explen - 2);
  auto rhsSig = (rhsConv.m_data.IEEE.significand | (one << siglen)) << (explen - 2);

  // calculate the difference between exponents
  int expDiff = lhsConv.m_data.IEEE.exponent - rhsConv.m_data.IEEE.exponent;
  auto& min = (expDiff > 0) ? rhsSig : lhsSig;
  expDiff *= signMap[!(expDiff > 0)];

  // get both operands in the same format with same exponent
  result.m_data.IEEE.exponent = std::max(lhsConv.m_data.IEEE.exponent, rhsConv.m_data.IEEE.exponent);
  min = roundingShift(min, expDiff);

  // use two's complement for negative
  lhsSig *= signMap[lhsConv.m_data.IEEE.sign];
  rhsSig *= signMap[rhsConv.m_data.IEEE.sign];

  // sum the significands for our result
  auto resultSig = lhsSig + rhsSig;

  // set the sign based on result
  // use two's complement on the result significand
  result.m_data.IEEE.sign = (resultSig >> (sizeof(TLargest) * CHAR_BIT - 1));
  resultSig *= signMap[result.m_data.IEEE.sign];

  // amount of overflow bits we need to normalize
  // this has to account for our previous bias
  auto overflow = msb(resultSig) - siglen - (explen - 1);

  // add the overflow to our exponent to normalize the float
  result.m_data.IEEE.exponent += overflow;
  // assign our calculated significand, if the max exponent has been reached we return infinity
  result.m_data.IEEE.significand = normalize<siglen>(resultSig) * (result.m_data.IEEE.exponent != maxExponent(explen));
  return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
template<unsigned TInSigLen, unsigned TInExpLen>
auto dynamicFloat<TSignificand, TExponent>::operator- (const dynamicFloat<TInSigLen, TInExpLen> _rhs) const
{
  return *this + (-_rhs);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
template<unsigned TInSigLen, unsigned TInExpLen>
dynamicFloat<TSignificand, TExponent>& dynamicFloat<TSignificand, TExponent>::operator += (const dynamicFloat<TInSigLen, TInExpLen> _rhs)
{
  *this = (*this) + _rhs;
  return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
template<typename F, typename>
dynamicFloat<TSignificand, TExponent>& dynamicFloat<TSignificand, TExponent>::operator += (const F _rhs)
{
  (*this) += detail::dynamicEquivalent<F>(_rhs);
  return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
template<unsigned TInSigLen, unsigned TInExpLen>
dynamicFloat<TSignificand, TExponent>& dynamicFloat<TSignificand, TExponent>::operator -= (const dynamicFloat<TInSigLen, TInExpLen> _rhs)
{
  *this = (*this) - _rhs;
  return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
template<typename F, typename>
dynamicFloat<TSignificand, TExponent>& dynamicFloat<TSignificand, TExponent>::operator -= (const F _rhs)
{
  (*this) -= detail::dynamicEquivalent<F>(_rhs);
  return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
template<unsigned TInSigLen, unsigned TInExpLen>
dynamicFloat<TSignificand, TExponent>& dynamicFloat<TSignificand, TExponent>::operator *= (const dynamicFloat<TInSigLen, TInExpLen> _rhs)
{
  *this = (*this) * _rhs;
  return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
template<typename F, typename>
dynamicFloat<TSignificand, TExponent>& dynamicFloat<TSignificand, TExponent>::operator *= (const F _rhs)
{
  (*this) *= detail::dynamicEquivalent<F>(_rhs);
  return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
template<unsigned TInSigLen, unsigned TInExpLen>
dynamicFloat<TSignificand, TExponent>& dynamicFloat<TSignificand, TExponent>::operator /= (const dynamicFloat<TInSigLen, TInExpLen> _rhs)
{
  *this = (*this) / _rhs;
  return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
template<typename F, typename>
dynamicFloat<TSignificand, TExponent>& dynamicFloat<TSignificand, TExponent>::operator /= (const F _rhs)
{
  (*this) /= detail::dynamicEquivalent<F>(_rhs);
  return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
dynamicFloat<TSignificand, TExponent> dynamicFloat<TSignificand, TExponent>::operator++(const int)
{
  auto copy = *this;
  (*this) += 1.0f;
  return copy;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
dynamicFloat<TSignificand, TExponent>& dynamicFloat<TSignificand, TExponent>::operator++()
{
  return (*this) += 1.0f;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
dynamicFloat<TSignificand, TExponent> dynamicFloat<TSignificand, TExponent>::operator--(const int)
{
  auto copy = *this;
  (*this) -= 1.0f;
  return copy;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
dynamicFloat<TSignificand, TExponent>& dynamicFloat<TSignificand, TExponent>::operator--()
{
  return (*this) -= 1.0f;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
dynamicFloat<TSignificand, TExponent> dynamicFloat<TSignificand, TExponent>::operator-() const
{
  dynamicFloat<TSignificand, TExponent> neg = *this;
  neg.m_data.IEEE.sign = ~neg.m_data.IEEE.sign;
  return neg;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//===============================================================================================================
//-------------------------Free functions for compatability with built in types----------------------------------
//===============================================================================================================

//----------------------------------------------Equality operator------------------------------------------------

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<unsigned TInSigLen, unsigned TInExpLen, typename F, typename = std::enable_if_t<std::is_floating_point<F>::value>>
bool operator== (const dynamicFloat<TInSigLen, TInExpLen> _lhs, const F _rhs)
{
  return _lhs == dynamicFloat<TInSigLen, TInExpLen>(_rhs);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<unsigned TInSigLen, unsigned TInExpLen, typename F, typename = std::enable_if_t<std::is_floating_point<F>::value>>
bool operator== (const F _lhs, const dynamicFloat<TInSigLen, TInExpLen> _rhs)
{
  return _rhs == dynamicFloat<TInSigLen, TInExpLen>(_lhs);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//-----------------------------------------------Addition operator-----------------------------------------------

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<unsigned TInSigLen, unsigned TInExpLen, typename F, typename = std::enable_if_t<std::is_floating_point<F>::value>>
auto operator+ (const dynamicFloat<TInSigLen, TInExpLen> _lhs, const F _rhs)
{
  return _lhs + detail::dynamicEquivalent<F>(_rhs);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<unsigned TInSigLen, unsigned TInExpLen, typename F, typename = std::enable_if_t<std::is_floating_point<F>::value>>
auto operator+ (const F _lhs, const dynamicFloat<TInSigLen, TInExpLen> _rhs)
{
  return _rhs + detail::dynamicEquivalent<F>(_lhs);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------Subtraction operator---------------------------------------------

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<unsigned TInSigLen, unsigned TInExpLen, typename F, typename = std::enable_if_t<std::is_floating_point<F>::value>>
auto operator- (const dynamicFloat<TInSigLen, TInExpLen> _lhs, const F _rhs)
{
  return _lhs - detail::dynamicEquivalent<F>(_rhs);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<unsigned TInSigLen, unsigned TInExpLen, typename F, typename = std::enable_if_t<std::is_floating_point<F>::value>>
auto operator- (const F _lhs, const dynamicFloat<TInSigLen, TInExpLen> _rhs)
{
  return _rhs - detail::dynamicEquivalent<F>(_lhs);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//--------------------------------------------Multiplication operator--------------------------------------------

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<unsigned TInSigLen, unsigned TInExpLen, typename F, typename = std::enable_if_t<std::is_floating_point<F>::value>>
auto operator* (const dynamicFloat<TInSigLen, TInExpLen> _lhs, const F _rhs)
{
  return _lhs * detail::dynamicEquivalent<F>(_rhs);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<unsigned TInSigLen, unsigned TInExpLen, typename F, typename = std::enable_if_t<std::is_floating_point<F>::value>>
auto operator* (const F _lhs, const dynamicFloat<TInSigLen, TInExpLen> _rhs)
{
  return _rhs * detail::dynamicEquivalent<F>(_lhs);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//-----------------------------------------------Division operator----------------------------------------------

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<unsigned TInSigLen, unsigned TInExpLen, typename F, typename = std::enable_if_t<std::is_floating_point<F>::value>>
auto operator/ (const dynamicFloat<TInSigLen, TInExpLen> _lhs, const F _rhs)
{
  return _lhs / detail::dynamicEquivalent<F>(_rhs);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<unsigned TInSigLen, unsigned TInExpLen, typename F, typename = std::enable_if_t<std::is_floating_point<F>::value>>
auto operator/ (const F _lhs, const dynamicFloat<TInSigLen, TInExpLen> _rhs)
{
  return _rhs / detail::dynamicEquivalent<F>(_lhs);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



//===============================================================================================================
//----------------------------------------------PRIVATE API------------------------------------------------------
//===============================================================================================================

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
bool dynamicFloat<TSignificand, TExponent>::isNaN() const noexcept
{
  return m_data.IEEE.significand != 0 && m_data.IEEE.exponent == maxExponent(TExponent);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
bool dynamicFloat<TSignificand, TExponent>::isInfinity() const noexcept
{
  return m_data.IEEE.significand == 0 && m_data.IEEE.exponent == maxExponent(TExponent);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
bool dynamicFloat<TSignificand, TExponent>::isDenorm() const noexcept
{
  return m_data.IEEE.exponent == 0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
constexpr auto dynamicFloat<TSignificand, TExponent>::bias(const unsigned _exponent) noexcept
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
constexpr auto dynamicFloat<TSignificand, TExponent>::maxExponent(const unsigned _exponent) noexcept
{
  return (1ul << _exponent) - 1;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
constexpr auto dynamicFloat<TSignificand, TExponent>::iShiftPow2(const unsigned _mantissa) noexcept
{
  return 1ul << _mantissa;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
constexpr unsigned dynamicFloat<TSignificand, TExponent>::msb(uint64_t v)
{
  constexpr int deBruijnTable[128] ={
     0, 48, -1, -1, 31, -1, 15, 51, -1, 63,  5, -1, -1, -1, 19, -1,
    23, 28, -1, -1, -1, 40, 36, 46, -1, 13, -1, -1, -1, 34, -1, 58,
    -1, 60,  2, 43, 55, -1, -1, -1, 50, 62,  4, -1, 18, 27, -1, 39,
    45, -1, -1, 33, 57, -1,  1, 54, -1, 49, -1, 17, -1, -1, 32, -1,
    53, -1, 16, -1, -1, 52, -1, -1, -1, 64,  6,  7,  8, -1,  9, -1,
    -1, -1, 20, 10, -1, -1, 24, -1, 29, -1, -1, 21, -1, 11, -1, -1,
    41, -1, 25, 37, -1, 47, -1, 30, 14, -1, -1, -1, -1, 22, -1, -1,
    35, 12, -1, -1, -1, 59, 42, -1, -1, 61,  3, 26, 38, 44, -1, 56
  };

  // is constant time
  for(int i = 1; i < 64; i <<= 1) v |= v >> i;

  return static_cast<unsigned>(deBruijnTable[static_cast<uint64_t>( v * 0x6c04f118e9966f6bUL ) >> 57]);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
constexpr unsigned dynamicFloat<TSignificand, TExponent>::lsb(const uint64_t v)
{
  constexpr int deBruijnTable[64] = {
    0,   1, 48,  2, 57, 49, 28,  3, 61, 58, 50, 42, 38, 29, 17,  4,
    62, 55, 59, 36, 53, 51, 43, 22, 45, 39, 33, 30, 24, 18, 12,  5,
    63, 47, 56, 27, 60, 41, 37, 16, 54, 35, 52, 21, 44, 32, 23, 11,
    46, 26, 40, 15, 34, 20, 31, 10, 25, 14, 19,  9, 13,  8,  7,  6
  };
  return static_cast<unsigned>(deBruijnTable[static_cast<uint64_t>((v & -v) * 0x03F79D71B4CB0A89U) >> 58]);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
template<int width, typename Tuint>
bool dynamicFloat<TSignificand, TExponent>::nonOverflowMult(const Tuint _lhs, const Tuint _rhs, Tuint &io_result) const
{
  // we need to calculate the product of all combnations between the high bits,
  // and the low bits of both operands. This will use long division to carry.

  constexpr unsigned size = sizeof(Tuint) * CHAR_BIT;
  constexpr unsigned halfSize = size / 2;
  // lambdas to get the lower and higher bits of an integer
  const auto lowBits = [](Tuint _in) { return ((Tuint{1} << halfSize) - 1) & _in; };
  const auto highBits = [](Tuint _in) { return _in >> halfSize; };

  // temp variable used to store intermediate results
  // calculate the low product
  auto temp = lowBits(_lhs) * lowBits(_rhs);
  // get the least signficant bits from the low product
  const auto lowestBits = lowBits(temp);

  // carry the high bits forward from the past calculation
  temp = highBits(_lhs) * lowBits(_rhs) + highBits(temp);
  // get the low, and the high to carry forward
  auto lowerMidBits = lowBits(temp);
  auto upperMidBits = highBits(temp);

  // add the next product to our low bits from last calc
  temp = lowerMidBits + lowBits(_lhs) * highBits(_rhs);
  // save the new low bits
  lowerMidBits = lowBits(temp);

  // repeat the previous calc using the new highest bits
  temp = upperMidBits + highBits(_lhs) * highBits(_rhs) + highBits(temp);
  upperMidBits = lowBits(temp);

  // get the final highest bits
  const auto highestBits = highBits(temp);

  // get the result produced from normal overflowing multiplication
  const auto overflowResult = lowerMidBits << halfSize | lowestBits;
  // calculate the bits that would be lost to overflow
  const auto carry = highestBits << halfSize | upperMidBits;

  // shift the overflow result to make room for carry
  // shift the carry to the front then add
  io_result = roundingShift(overflowResult, msb(carry)) + (carry << (size - msb(carry)));

  // is the carry empty?
  const bool emptyCarry = !!carry;

  // calculate the number of digits that the result occupies
  // if the carry isn't empty we need to add 64 as it is at least this long
  // otherwise we get the number of digits from the first word
  const unsigned numDigits = msb(carry) + !emptyCarry * msb(overflowResult) + emptyCarry * size;

  // return whether the exponent requires an increment
  return numDigits >= width;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
template <int width, typename Tuint>
bool dynamicFloat<TSignificand, TExponent>::longDivide(Tuint _lhs, const Tuint _rhs, Tuint &io_result) const
{
  constexpr unsigned size = sizeof(Tuint) * CHAR_BIT;
  Tuint one = 1;

  unsigned shift = size - msb(_lhs);
  // shift the left operand left the max amount
  // to get as much precision as posible
  _lhs <<= shift;
  bool set = false;
  Tuint rem = io_result = 0;
  // start from the second lest significant bit
  int i = size-1;
  for(int p = i; p >= 0; --i, rem <<= 1)
  {
    // we reverse the shift i is negative
    bool neg = i < 0;
    auto pos = ((one >> -i * neg) << i * !neg);
    // get the bit stored at the current bit
    rem |= !!(_lhs & pos);
    // if the remainder is larger than the right operand
    // we subtract and store the bit
    bool subtract = (rem >= _rhs);
    io_result |= (one << p) * subtract;
    rem -= _rhs * subtract;
    // marks that we have started storing bits
    set |= subtract;
    // only decrement the postion to store bits
    // after the first has been stored, so that
    // max precision is used, and order preserved
    p -= set;
  }
  // normalize the result
  io_result = roundingShift(io_result, shift);
  return i < -width;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
template <typename Tuint>
auto dynamicFloat<TSignificand, TExponent>::roundingShift(Tuint io_num, unsigned _shift) const
{
  // check whether the shift exceeds the type width
  unsigned notExceeded = !(_shift >= sizeof(Tuint) * CHAR_BIT);
  Tuint one = 1;
  // save the bits that will be lost by shifting
  Tuint loss = io_num & (((one << _shift) - 1) * notExceeded);
  // check whether the shift is negative
  bool neg = _shift >> (sizeof(_shift) * CHAR_BIT - 1);
  // shift in the opposite direction for negative
  io_num >>= (_shift * !neg);
  io_num <<= (-_shift * neg);
  // we multiply to set to zero if the shift was same width as type,
  // this is required as shifting by the width of the type is UB,
  // but we would expect a value of zero as the result.
  io_num *= notExceeded;
  // avoid branch by adding the bool
  io_num += static_cast<unsigned>(
        _shift && notExceeded && !neg &&                       // check that the shift is non-zero and not too large
        (msb(loss) == _shift)                                  // check that first bit is 1
        && ((loss & ~(one << (_shift - 1))) || (io_num & 1) )  // check that there are other 1's
        );
  return io_num;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
template <unsigned width, typename Tuint>
auto dynamicFloat<TSignificand, TExponent>::normalize(Tuint io_num) const
{
  // amount of overflow bits we need to normalize
  auto overflow = msb(io_num) - width - 1;
  // check sign of new significand
  if(overflow >> (sizeof(overflow) * CHAR_BIT - 1))
  {
    // don't need to round, reverse shift for negative
    io_num <<= -overflow;
  }
  // check not zero
  else
  {
    io_num = roundingShift(io_num, overflow);
  }
  return io_num;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#endif
