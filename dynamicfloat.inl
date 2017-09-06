
#ifndef DYNAMICFLOAT_INL_INCLUDED
#define DYNAMICFLOAT_INL_INCLUDED


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
auto dynamicFloat<TSignificand, TExponent>::operator* (const dynamicFloat<TInSigLen, TInExpLen> _rhs) const
{
  using TLargest = typename std::conditional_t< ((TSignificand + TExponent) > (TInSigLen + TInExpLen)), dynamicFloat<TSignificand, TExponent>, dynamicFloat<TInSigLen, TInExpLen>>;
  constexpr unsigned siglen = std::max(TSignificand,TInSigLen);
  constexpr unsigned explen = std::max(TExponent,TInExpLen);
  TLargest lhsConv = *this;
  TLargest rhsConv = _rhs;

  //if either is zero return the other
  if(!lhsConv.m_data.m_bits || !rhsConv.m_data.m_bits) return TLargest{0.f};

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
  // m_data.IEEE uses 1.xxx... where the 1 isn't stored, so we add it
  auto lhsSig = lhsConv.m_data.IEEE.significand | (1ul << siglen);
  auto rhsSig = rhsConv.m_data.IEEE.significand | (1ul << siglen);

  // exponent is the sum of the operands exponents
  result.m_data.IEEE.exponent = (lhsConv.m_data.IEEE.exponent - bias(explen)) + (rhsConv.m_data.IEEE.exponent - bias(explen)) + bias(explen);

  uint64_t resultSig;
  bool up = nonOverflowMult<(siglen + 1) * 2>(lhsSig, rhsSig, resultSig);

  // amount of overflow bits we need to normalize
  auto overflow = msb(resultSig) - siglen - 1;

  // check sign of new significand
  if(overflow >> (sizeof(overflow) * 8 - 1))
  {
    // don't need to round, reverse shift for negative
    resultSig <<= -overflow;
  }
  // check not zero
  else if (overflow)
  {
    resultSig = roundingShift(resultSig, overflow);
  }

  // XOR sign
  result.m_data.IEEE.sign = (lhsConv.m_data.IEEE.sign != rhsConv.m_data.IEEE.sign);
  // assign our calculated significand
  result.m_data.IEEE.significand = resultSig;
  // avoid a branch here and just add the bool
  result.m_data.IEEE.exponent += static_cast<int>(up);

  return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
template<unsigned TInSigLen, unsigned TInExpLen>
auto dynamicFloat<TSignificand, TExponent>::operator/ (const dynamicFloat<TInSigLen, TInExpLen> _rhs) const
{
  using TLargest = typename std::conditional_t< ((TSignificand + TExponent) > (TInSigLen + TInExpLen)), dynamicFloat<TSignificand, TExponent>, dynamicFloat<TInSigLen, TInExpLen>>;
  constexpr unsigned siglen = std::max(TSignificand,TInSigLen);
  constexpr unsigned explen = std::max(TExponent,TInExpLen);
  TLargest lhsConv = *this;
  TLargest rhsConv = _rhs;

  //if either is zero return the other
  if(!lhsConv.m_data.m_bits || !rhsConv.m_data.m_bits) return TLargest{0.f};

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
  // m_data.IEEE uses 1.xxx... where the 1 isn't stored, so we add it
  auto lhsSig = lhsConv.m_data.IEEE.significand | (1ul << siglen);
  auto rhsSig = rhsConv.m_data.IEEE.significand | (1ul << siglen);

  // exponent is the sum of the operands exponents
  result.m_data.IEEE.exponent = (lhsConv.m_data.IEEE.exponent - bias(explen)) - (rhsConv.m_data.IEEE.exponent - bias(explen)) + bias(explen);
  // sum the significands for our result
  uint64_t resultSig;

  bool down = longDivide<siglen + 1>(lhsSig , rhsSig, resultSig);

  // amount of overflow bits we need to normalize
  auto overflow = msb(resultSig) - siglen - 1;
  std::cout<<overflow<<'\n';

  // check sign of new significand
  if(overflow >> (sizeof(overflow) * 8 - 1))
  {
    // don't need to round, reverse shift for negative
    resultSig <<= -overflow;
  }
  // check not zero
  else if (overflow)
  {
    resultSig = roundingShift(resultSig, overflow);
  }

  // XOR sign
  result.m_data.IEEE.sign = (lhsConv.m_data.IEEE.sign != rhsConv.m_data.IEEE.sign);
  // assign our calculated significand
  result.m_data.IEEE.significand = resultSig;

  result.m_data.IEEE.exponent -= static_cast<int>(down);

  return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
template<unsigned TInSigLen, unsigned TInExpLen>
auto dynamicFloat<TSignificand, TExponent>::operator+ (const dynamicFloat<TInSigLen, TInExpLen> _rhs) const
{
  // define the return type as the most precision passed
  using TLargest = typename std::conditional_t< ((TSignificand + TExponent) > (TInSigLen + TInExpLen)), dynamicFloat<TSignificand, TExponent>, dynamicFloat<TInSigLen, TInExpLen>>;

  // if either passed is NaN the result is NaN
  if(isNaN() || _rhs.isNaN())
  {
    std::cout<<"CAUGHT\n";
    return std::numeric_limits<TLargest>::signaling_NaN();
  }

  // if either passed is infinity the result is infinity
  else if(isInfinity() || _rhs.isInfinity())
  {
    return std::numeric_limits<TLargest>::infinity();
  }

  constexpr unsigned siglen = std::max(TSignificand,TInSigLen);
  constexpr unsigned explen = std::max(TExponent,TInExpLen);
  TLargest lhsConv = *this;
  TLargest rhsConv = _rhs;

  //if either is zero return the other
  if(!lhsConv.m_data.m_bits) return rhsConv;
  if(!rhsConv.m_data.m_bits) return lhsConv;

  TLargest result;
  // calculate the difference between exponents
  int expDiff = lhsConv.m_data.IEEE.exponent - rhsConv.m_data.IEEE.exponent;

  // get the significands and add the hidden bit,
  // m_data.IEEE uses 1.xxx... where the 1 isn't stored, so we add it
  auto lhsSig = lhsConv.m_data.IEEE.significand | (1ul << siglen);
  auto rhsSig = rhsConv.m_data.IEEE.significand | (1ul << siglen);

  // get both operands in the same format with same exponent
  if(expDiff > 0)
  {
    // choose this format's exponent
    result.m_data.IEEE.exponent = lhsConv.m_data.IEEE.exponent;
    rhsSig = roundingShift(rhsSig, expDiff);
  }
  else
  {
    expDiff = -expDiff;
    // choose this format's exponent
    result.m_data.IEEE.exponent = rhsConv.m_data.IEEE.exponent;
    lhsSig = roundingShift(lhsSig, expDiff);
  }

  // use two's complement for negative
  if(lhsConv.m_data.IEEE.sign) lhsSig = -lhsSig;
  if(rhsConv.m_data.IEEE.sign) rhsSig = -rhsSig;

  // sum the significands for our result
  auto resultSig = lhsSig + rhsSig;

  // path dependent on sign of result
  if(resultSig >> (sizeof(TLargest) * 8 - 1))
  {
    // negative sign
    result.m_data.IEEE.sign = 1;
    // convert result significand to two's complement
    result.m_data.IEEE.significand = -resultSig;
  }
  else
  {
    // positive sign
    result.m_data.IEEE.sign = 0;
    // amount of overflow bits we need to normalize
    auto overflow = msb(resultSig) - siglen - 1;

    // check sign of new significand
    if(overflow >> (sizeof(overflow) * 8 - 1))
    {
      // don't need to round, reverse shift for negative
      resultSig <<= -overflow;
    }
    // check not zero
    else if (overflow)
    {
      resultSig = roundingShift(resultSig, overflow);
    }

    // add the overflow to our exponent to normalize the float
    result.m_data.IEEE.exponent += overflow;
    if(result.m_data.IEEE.exponent != maxExponent(explen))
      // assign our calculated significand
      result.m_data.IEEE.significand = resultSig;
    else
      result.m_data.IEEE.significand = 0;

  }
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
dynamicFloat<TSignificand, TExponent> dynamicFloat<TSignificand, TExponent>::operator--(const int) const
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
  neg.m_data.IEEE.sign = ~neg.m_data.IEEE.sign;
  return neg;
}
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

  return static_cast<unsigned>(deBruijnTable[static_cast<uint64_t>( v * 0x6c04f118e9966f6bUL ) >> 57]);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
constexpr unsigned dynamicFloat<TSignificand, TExponent>::lsb(const uint64_t v)
{
  constexpr int deBruijnTable[64] =
  {
    0, 1, 48, 2, 57, 49, 28, 3, 61, 58, 50, 42, 38, 29, 17, 4,
    62, 55, 59, 36, 53, 51, 43, 22, 45, 39, 33, 30, 24, 18, 12, 5,
    63, 47, 56, 27, 60, 41, 37, 16, 54, 35, 52, 21, 44, 32, 23, 11,
    46, 26, 40, 15, 34, 20, 31, 10, 25, 14, 19, 9, 13, 8, 7, 6
  };
  return static_cast<unsigned>(deBruijnTable[static_cast<uint64_t>((v & -v) * 0x03F79D71B4CB0A89U) >> 58]);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
template<int width>
bool dynamicFloat<TSignificand, TExponent>::nonOverflowMult(const uint64_t _lhs, const uint64_t _rhs, uint64_t &io_result) const
{
  // we need to calculate the product of all combnations between the high bits,
  // and the low bits of both operands. This will use long division to carry.

  // lambdas to get the lower and higher bits of an integer
  const auto lowBits = [](uint64_t _in) { return ((1ul << 32) - 1) & _in; };
  const auto highBits = [](uint64_t _in) { return _in >> 32; };

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
  const auto overflowResult = lowerMidBits << 32 | lowestBits;
  // calculate the bits that would be lost to overflow
  const auto carry = highestBits << 32 | upperMidBits;

  // shift the overflow result to make room for carry
  // shift the carry to the front then add
  io_result = roundingShift(overflowResult, msb(carry)) + (carry << (64 - msb(carry)));

  // is the carry empty?
  const bool emptyCarry = static_cast<bool>(carry);

  // calculate the number of digits that the result occupies
  // if the carry isn't empty we need to add 64 as it is at least this long
  // otherwise we get the number of digits from the first word
  const unsigned numDigits = msb(carry) + static_cast<unsigned>(!emptyCarry) * msb(overflowResult) + static_cast<unsigned>(emptyCarry) * 64;

  // return whether the exponent requires an increment
  return numDigits >= width;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
template <int width>
bool dynamicFloat<TSignificand, TExponent>::longDivide(uint64_t _lhs, const uint64_t _rhs, uint64_t &io_result) const
{
  unsigned shift = 64 - msb(_lhs);
  _lhs <<= shift;
  bool set = false;
  io_result = 0;
  uint64_t rem = 0;
  int i = 63;
  for(int p = 63; p >= 0; --i)
  {
    rem <<= 1;
    if(i >= 0)
      rem |= static_cast<unsigned>(static_cast<bool>(_lhs & (1ul << i)));
    else
      rem |= static_cast<unsigned>(static_cast<bool>(_lhs & (1ul >> -i)));
    if(rem >= _rhs)
    {
      rem -= _rhs;
      io_result |= 1ul << p;
      set = true;
    }
    p -= static_cast<int>(set);
  }
  io_result = roundingShift(io_result, shift);

  return i < -width;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned TSignificand, unsigned TExponent>
auto dynamicFloat<TSignificand, TExponent>::roundingShift(uint64_t io_num, const uint64_t _shift) const
{
  if(!_shift) return io_num;
  // save the bits that will be lost by shifting
  auto loss = io_num & ((1ul << _shift) - 1);
  io_num >>= _shift;
  if( (msb(loss) == static_cast<unsigned>(_shift))           // check that first bit is 1
      && ((loss & ~(1ul << (_shift - 1))) || (io_num & 1) )) // check that there are other 1's
  {
    io_num++;
  }
  return io_num;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#endif
