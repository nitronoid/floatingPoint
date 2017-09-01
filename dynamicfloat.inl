
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
  return *this;
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
auto dynamicFloat<TSignificand, TExponent>::operator* (const dynamicFloat<TInSigLen, TInExpLen> _rhs)
{
  using TLargest = typename std::conditional_t< ((TSignificand + TExponent) > (TInSigLen + TInExpLen)), dynamicFloat<TSignificand, TExponent>, dynamicFloat<TInSigLen, TInExpLen>>;
  constexpr unsigned siglen = std::max(TSignificand,TInSigLen);
  constexpr unsigned explen = std::max(TExponent,TInExpLen);
  TLargest lhsConv = *this;
  TLargest rhsConv { _rhs};

  //if either is zero return the other
  if(!lhsConv.m_bits || !rhsConv.m_bits) return TLargest{0.f};

  // if either passed is NaN the result is NaN
  if(lhsConv.isNaN() || rhsConv.isNaN())
  {
    return std::numeric_limits<TLargest>::signaling_NaN();
  }

  // if either passed is infinity the result is infinity
  else if(lhsConv.isInfinity() || rhsConv.isInfinity())
  {
    return std::numeric_limits<TLargest>::infinity();
  }

  TLargest result;

  // get the significands and add the hidden bit,
  // IEEE uses 1.xxx... where the 1 isn't stored, so we add it
  auto lhsSig = lhsConv.IEEE.significand | (uint64_t{1} << siglen);
  auto rhsSig = rhsConv.IEEE.significand | (uint64_t{1} << siglen);

  // remove trailing zeros to prevent overflow durng multiplcation
  // automatically corrected when we normalize
  lhsSig >>= lsb(lhsSig);
  rhsSig >>= lsb(rhsSig);


  // exponent is the sum of the operands exponents
  result.IEEE.exponent = (rhsConv.IEEE.exponent-bias(explen)) + (lhsConv.IEEE.exponent-bias(explen)) + bias(explen);

  // sum the significands for our result
  auto resultSig = lhsSig * rhsSig;

  // amount of overflow bits we need to normalize
  auto overflow = msb(resultSig) - siglen;

  //calculate the required exponent change after normalize
  auto normExpChange = msb(resultSig) - (msb(lhsSig) + msb(rhsSig));
  std::cout<<normExpChange<<'\n';
  // check sign of new significand
  if(overflow >> (sizeof(overflow) * 8 - 1))
  {
    // don't need to round, reverse shift for negative
    resultSig <<= -overflow;
  }
  // check not zero
  else if (overflow)
  {
    // the bits lost when the significand is normalized
//    auto loss = resultSig & ((uint64_t{1} << overflow) - 1);
    // shift with rounding
    std::cout<<std::bitset<32>((resultSig & (uint64_t{1}<<(overflow-1))))<<'\n';
    resultSig = (resultSig >> overflow) + msb((resultSig & (uint64_t{1}<<(overflow-1))) << 1);

//    // round the signficand
//    if( (msb(loss) == static_cast<unsigned>(overflow - 1)) // check that first bit is 1
//        && (loss & ~(1<<(overflow-1))))                     // check that there are other 1's
//    {
//      resultSig++;
//    }
  }




  // XOR sign
  result.IEEE.sign = (lhsConv.IEEE.sign != rhsConv.IEEE.sign);
  // assign our calculated significand
  result.IEEE.significand = resultSig;
  // correct the exponent
  result.IEEE.exponent += normExpChange;
  return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
template<unsigned TInSigLen, unsigned TInExpLen>
auto dynamicFloat<TSignificand, TExponent>::operator+ (const dynamicFloat<TInSigLen, TInExpLen> _rhs)
{
  // define the return type as the most precision passed
  using TLargest = typename std::conditional_t< ((TSignificand + TExponent) > (TInSigLen + TInExpLen)), dynamicFloat<TSignificand, TExponent>, dynamicFloat<TInSigLen, TInExpLen>>;
  constexpr unsigned siglen = std::max(TSignificand,TInSigLen);
  TLargest lhsConv = *this;
  TLargest rhsConv { _rhs};

  //if either is zero return the other
  if(!lhsConv.m_bits) return rhsConv;
  if(!rhsConv.m_bits) return lhsConv;

  // if either passed is NaN the result is NaN
  if(lhsConv.isNaN() || rhsConv.isNaN())
  {
    return std::numeric_limits<TLargest>::signaling_NaN();
  }

  // if either passed is infinity the result is infinity
  else if(lhsConv.isInfinity() || rhsConv.isInfinity())
  {
    return std::numeric_limits<TLargest>::infinity();
  }

  TLargest result;
  // calculate the difference between exponents
  int expDiff = lhsConv.IEEE.exponent - rhsConv.IEEE.exponent;

  // get the significands and add the hidden bit,
  // IEEE uses 1.xxx... where the 1 isn't stored, so we add it
  auto lhsSig = lhsConv.IEEE.significand | (uint64_t{1} << siglen);
  auto rhsSig = rhsConv.IEEE.significand | (uint64_t{1} << siglen);
  uint64_t loss;

  // get both operands in the same format with same exponent
  if(expDiff > 0)
  {
    loss = rhsSig & ((uint64_t{1} << expDiff) - 1);
    rhsSig >>= expDiff;
    // choose this format's exponent
    result.IEEE.exponent = lhsConv.IEEE.exponent;
  }
  else
  {
    expDiff = -expDiff;
    loss = lhsSig & ((uint64_t{1} << expDiff) - 1);
    lhsSig >>= expDiff;
    // choose this format's exponent
    result.IEEE.exponent = rhsConv.IEEE.exponent;
  }

  if( (msb(loss) == static_cast<unsigned>(expDiff - 1)) // check that first bit is 1
      && (loss & ~(1<<(expDiff-1))))                     // check that there are other 1's
    rhsSig++;

  // use two's complement for negative
  if(lhsConv.IEEE.sign) lhsSig = -lhsSig;
  if(rhsConv.IEEE.sign) rhsSig = -rhsSig;

  // sum the significands for our result
  auto resultSig = lhsSig + rhsSig;

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
      // don't need to round, reverse shift for negative
      resultSig <<= -overflow;
    // check not zero
    else if (overflow)
      // shift with rounding
      resultSig = (resultSig >> overflow) + msb((resultSig & (uint64_t{1}<<(overflow-1))) << 1);

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
template<unsigned TInSigLen, unsigned TInExpLen>
auto dynamicFloat<TSignificand, TExponent>::operator- (const dynamicFloat<TInSigLen, TInExpLen> _rhs)
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
template<unsigned TInSigLen, unsigned TInExpLen>
dynamicFloat<TSignificand, TExponent>& dynamicFloat<TSignificand, TExponent>::operator -= (const dynamicFloat<TInSigLen, TInExpLen> _rhs)
{
  *this = (*this) - _rhs;
  return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//template <unsigned TSignificand, unsigned TExponent>
//template<unsigned TInSigLen, unsigned TInExpLen>
//dynamicFloat<TSignificand, TExponent>& dynamicFloat<TSignificand, TExponent>::operator *= (const dynamicFloat<TInSigLen, TInExpLen> _rhs)
//{

//}
////-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//template <unsigned TSignificand, unsigned TExponent>
//template<unsigned TInSigLen, unsigned TInExpLen>
//dynamicFloat<TSignificand, TExponent>& dynamicFloat<TSignificand, TExponent>::operator /= (const dynamicFloat<TInSigLen, TInExpLen> _rhs)
//{

//}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
dynamicFloat<TSignificand, TExponent> dynamicFloat<TSignificand, TExponent>::operator--(int) const
{

}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
dynamicFloat<TSignificand, TExponent>& dynamicFloat<TSignificand, TExponent>::operator--() const
{

}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
dynamicFloat<TSignificand, TExponent> dynamicFloat<TSignificand, TExponent>::operator-() const
{
  dynamicFloat<TSignificand, TExponent> neg = *this;
  neg.IEEE.sign = ~neg.IEEE.sign;
  return neg;
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
constexpr unsigned dynamicFloat<TSignificand, TExponent>::msb(uint64_t v)
{
  constexpr int deBruijnTable[128] ={
    0, // change to 1 if you want bitSize(0) = 1
    48, -1, -1, 31, -1, 15, 51, -1, 63, 5, -1, -1, -1, 19, -1,
    23, 28, -1, -1, -1, 40, 36, 46, -1, 13, -1, -1, -1, 34, -1, 58,
    -1, 60, 2, 43, 55, -1, -1, -1, 50, 62, 4, -1, 18, 27, -1, 39,
    45, -1, -1, 33, 57, -1, 1, 54, -1, 49, -1, 17, -1, -1, 32, -1,
    53, -1, 16, -1, -1, 52, -1, -1, -1, 64, 6, 7, 8, -1, 9, -1,
    -1, -1, 20, 10, -1, -1, 24, -1, 29, -1, -1, 21, -1, 11, -1, -1,
    41, -1, 25, 37, -1, 47, -1, 30, 14, -1, -1, -1, -1, 22, -1, -1,
    35, 12, -1, -1, -1, 59, 42, -1, -1, 61, 3, 26, 38, 44, -1, 56
  };

  // is constant time
  for(int i = 1; i < 64; i <<= 1) v |= v >> i;

  return deBruijnTable[static_cast<uint64_t>( v * 0x6c04f118e9966f6bUL ) >> 57] - 1;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
constexpr unsigned dynamicFloat<TSignificand, TExponent>::lsb(uint64_t v)
{
  constexpr int deBruijnTable[32] =
  {
    0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
    31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
  };
  return deBruijnTable[((uint32_t)((v & -v) * 0x077CB531U)) >> 27];
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#endif
