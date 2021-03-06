#include <cstdlib>
#include <iostream>
#include <cmath>
#include <bitset>
#include "dynamicfloat.h"
#include <QString>
#include <QtTest>
#include <random>

#define USE_DOUBLE

using float8 = dynamicFloat<3,4>;
using float16 = dynamicFloat<10,5>;
using float24 = dynamicFloat<16,7>;
using float32 = dynamicFloat<23,8>;
using float40 = dynamicFloat<31,8>;
using float48 = dynamicFloat<39,8>;
using float64 = dynamicFloat<52,11>;

class QDynamicFloatTest : public QObject
{
  Q_OBJECT

#ifdef USE_DOUBLE
  using TBuiltin = double;
  using TDynaFloat = float64;
  using Tuint = uint64_t;
  static constexpr unsigned width = 64;
#else
  using TBuiltin = float;
  using TDynaFloat = float32;
  using Tuint = uint32_t;
  static constexpr unsigned width = 32;
#endif

  static constexpr int iter = 100000;

public:
  QDynamicFloatTest() = default;

  std::random_device rd;
  std::mt19937 e2{rd()};
  std::uniform_real_distribution<TBuiltin> dist{static_cast<float>(-1000000000000000000), static_cast<float>(1000000000000000000)};

private Q_SLOTS:

  // Equality tests
  void testCaseEquality();

  // Arithmetic tests
  void testCaseAddition();
  void testCaseMultiplication();
  void testCaseDivision();
  void incrementDecrementTestCase();

  // Self arithmetic tests
  void selfAdditionTestCase();
  void selfAdditionBuiltinTestCase();
  void selfSubtractionTestCase();
  void selfSubtractionBuiltinTestCase();
  void selfMultiplicationTestCase();
  void selfMultiplicationBuiltinTestCase();
  void selfDivisionTestCase();
  void selfDivisionBuiltinTestCase();

  void foo();
};

void QDynamicFloatTest::foo()
{
  TBuiltin foo = 5.7771e+17;
  TBuiltin bar = -5.75865e+17;
  TBuiltin res = foo + bar;
  TDynaFloat food = foo;
  TDynaFloat bard = bar;
  TDynaFloat resd = food + bard;
  QCOMPARE(TBuiltin(resd), res);
}

void QDynamicFloatTest::testCaseEquality()
{
  TBuiltin builtinZero = 0.0f;
  TDynaFloat dynamicZero = 0.0f;
  TDynaFloat dynamicZeroConv = builtinZero;
  QVERIFY2(builtinZero == dynamicZero, "Failure");
  QVERIFY2(builtinZero == dynamicZeroConv, "Failure");

  TBuiltin builtinInf  = std::numeric_limits<TBuiltin>::infinity();
  TDynaFloat dynamicInf  = std::numeric_limits<TDynaFloat>::infinity();
  TDynaFloat dynamicInfConv = builtinInf;
  QVERIFY2(builtinInf == dynamicInf, "Failure");
  QVERIFY2(builtinInf == dynamicInfConv, "Failure");

  TBuiltin builtinSNaN = std::numeric_limits<TBuiltin>::signaling_NaN();
  TDynaFloat dynamicSNaN = std::numeric_limits<TDynaFloat>::signaling_NaN();
  TDynaFloat dynamicSNaNConv = builtinSNaN;
  QVERIFY2(
        !(dynamicSNaN == builtinSNaN) &&
        (std::bitset<width>(detail::bit_cast<Tuint>(builtinSNaN)) == std::bitset<width>(detail::bit_cast<Tuint>(dynamicSNaN))),
        "Failure"
        );
  QVERIFY2(
        !(dynamicSNaNConv == builtinSNaN) &&
        (std::bitset<width>(detail::bit_cast<Tuint>(builtinSNaN)) == std::bitset<width>(detail::bit_cast<Tuint>(dynamicSNaNConv))),
        "Failure"
        );

  TBuiltin builtinQNaN = std::numeric_limits<TBuiltin>::quiet_NaN();
  TDynaFloat dynamicQNaN = std::numeric_limits<TDynaFloat>::quiet_NaN();
  TDynaFloat dynamicQNaNConv = builtinQNaN;
  QVERIFY2(
        !(dynamicQNaN == builtinQNaN) &&
        (std::bitset<width>(detail::bit_cast<Tuint>(builtinQNaN)) == std::bitset<width>(detail::bit_cast<Tuint>(dynamicQNaN))),
        "Failure"
        );
  QVERIFY2(
        !(dynamicQNaNConv == builtinQNaN) &&
        (std::bitset<width>(detail::bit_cast<Tuint>(builtinQNaN)) == std::bitset<width>(detail::bit_cast<Tuint>(dynamicQNaNConv))),
        "Failure"
        );

  TBuiltin builtinMin  = std::numeric_limits<TBuiltin>::min();
  TDynaFloat dynamicMin  = std::numeric_limits<TDynaFloat>::min();
  TDynaFloat dynamicMinConv  = builtinMin;
  QVERIFY2(builtinMin == dynamicMin, "Failure");
  QVERIFY2(builtinMin == dynamicMinConv, "Failure");

  TBuiltin builtinMax  = std::numeric_limits<TBuiltin>::max();
  TDynaFloat dynamicMax  = std::numeric_limits<TDynaFloat>::max();
  TDynaFloat dynamicMaxConv  = builtinMax;
  QVERIFY2(builtinMax == dynamicMax, "Failure");
  QVERIFY2(builtinMax == dynamicMaxConv, "Failure");

  TBuiltin builtinLowest  = std::numeric_limits<TBuiltin>::lowest();
  TDynaFloat dynamicLowest  = std::numeric_limits<TDynaFloat>::lowest();
  TDynaFloat dynamicLowestConv  = builtinLowest;
  QVERIFY2(builtinLowest == dynamicLowest, "Failure");
  QVERIFY2(builtinLowest == dynamicLowestConv, "Failure");

  for(int i = 0; i < iter; ++i)
  {

    TBuiltin builtin = dist(e2);
    TDynaFloat dynamic = builtin;
    QVERIFY2(builtin == dynamic, "Failure");
  }
}

void QDynamicFloatTest::testCaseAddition()
{
  for(int i = 0; i < 100; ++i)
  {
    TBuiltin builtinLhs = dist(e2);
    TBuiltin builtinRhs = dist(e2);
    TBuiltin builtinResult = builtinLhs + builtinRhs;
    TDynaFloat dynamicLhs = builtinLhs;
    TDynaFloat dynamicRhs = builtinRhs;
    TDynaFloat dynamicResult = dynamicLhs + dynamicRhs;
    QCOMPARE(TBuiltin(dynamicResult), builtinResult);
  }
}

void QDynamicFloatTest::testCaseMultiplication()
{
  for(int i = 0; i < iter; ++i)
  {
    TBuiltin builtinLhs = dist(e2);
    TBuiltin builtinRhs = dist(e2);
    TBuiltin builtinResult = builtinLhs * builtinRhs;
    TDynaFloat dynamicLhs = builtinLhs;
    TDynaFloat dynamicRhs = builtinRhs;
    TDynaFloat dynamicResult = dynamicLhs * dynamicRhs;
    QCOMPARE(TBuiltin(dynamicResult), builtinResult);
  }
}

void QDynamicFloatTest::testCaseDivision()
{
  for(int i = 0; i < iter; ++i)
  {
    TBuiltin builtinLhs = dist(e2);
    TBuiltin builtinRhs = dist(e2);
    TBuiltin builtinResult = builtinLhs / builtinRhs;
    TDynaFloat dynamicLhs = builtinLhs;
    TDynaFloat dynamicRhs = builtinRhs;
    TDynaFloat dynamicResult = dynamicLhs / dynamicRhs;
    QCOMPARE(TBuiltin(dynamicResult), builtinResult);
  }
}

void QDynamicFloatTest::selfAdditionTestCase()
{
  for(int i = 0; i < iter; ++i)
  {
    //Addition
    TBuiltin builtinLhs = dist(e2);
    TBuiltin builtinRhs = dist(e2);

    TDynaFloat dynamicLhs = builtinLhs;
    TDynaFloat dynamicRhs = builtinRhs;
    builtinLhs += builtinRhs;
    dynamicLhs += dynamicRhs;
    QCOMPARE(TBuiltin(dynamicLhs), builtinLhs);
  }
}

void QDynamicFloatTest::selfAdditionBuiltinTestCase()
{
  for(int i = 0; i < iter; ++i)
  {
    //Addition with builtin
    TBuiltin builtinLhs = dist(e2);
    TBuiltin builtinRhs = dist(e2);

    TDynaFloat dynamicLhs = builtinLhs;

    builtinLhs += builtinRhs;
    dynamicLhs += builtinRhs;
    QCOMPARE(TBuiltin(dynamicLhs), builtinLhs);
  }
}

void QDynamicFloatTest::selfSubtractionTestCase()
{
  for(int i = 0; i < iter; ++i)
  {
    //Subtraction
    TBuiltin builtinLhs = dist(e2);
    TBuiltin builtinRhs = dist(e2);

    TDynaFloat dynamicLhs = builtinLhs;
    TDynaFloat dynamicRhs = builtinRhs;

    builtinLhs -= builtinRhs;
    dynamicLhs -= dynamicRhs;
    QCOMPARE(TBuiltin(dynamicLhs), builtinLhs);
  }
}

void QDynamicFloatTest::selfSubtractionBuiltinTestCase()
{
  for(int i = 0; i < iter; ++i)
  {
    //Subtraction with builtin
    TBuiltin builtinLhs = dist(e2);
    TBuiltin builtinRhs = dist(e2);

    TDynaFloat dynamicLhs = builtinLhs;

    builtinLhs -= builtinRhs;
    dynamicLhs -= builtinRhs;
    QCOMPARE(TBuiltin(dynamicLhs), builtinLhs);
  }
}

void QDynamicFloatTest::selfMultiplicationTestCase()
{
  for(int i = 0; i < iter; ++i)
  {
    //Multiplication
    TBuiltin builtinLhs = dist(e2);
    TBuiltin builtinRhs = dist(e2);

    TDynaFloat dynamicLhs = builtinLhs;
    TDynaFloat dynamicRhs = builtinRhs;

    builtinLhs *= builtinRhs;
    dynamicLhs *= dynamicRhs;
    QCOMPARE(TBuiltin(dynamicLhs), builtinLhs);
  }
}

void QDynamicFloatTest::selfMultiplicationBuiltinTestCase()
{
  for(int i = 0; i < iter; ++i)
  {
    //Multiplication with builtin
    TBuiltin builtinLhs = dist(e2);
    TBuiltin builtinRhs = dist(e2);

    TDynaFloat dynamicLhs = builtinLhs;

    builtinLhs *= builtinRhs;
    dynamicLhs *= builtinRhs;
    QCOMPARE(TBuiltin(dynamicLhs), builtinLhs);
  }
}

void QDynamicFloatTest::selfDivisionTestCase()
{
  for(int i = 0; i < iter; ++i)
  {
    //Division
    TBuiltin builtinLhs = dist(e2);
    TBuiltin builtinRhs = dist(e2);

    TDynaFloat dynamicLhs = builtinLhs;
    TDynaFloat dynamicRhs = builtinRhs;

    builtinLhs /= builtinRhs;
    dynamicLhs /= dynamicRhs;
    QCOMPARE(TBuiltin(dynamicLhs), builtinLhs);
  }
}

void QDynamicFloatTest::selfDivisionBuiltinTestCase()
{
  for(int i = 0; i < iter; ++i)
  {
    //Division with builtin
    TBuiltin builtinLhs = dist(e2);
    TBuiltin builtinRhs = dist(e2);

    TDynaFloat dynamicLhs = builtinLhs;

    builtinLhs /= builtinRhs;
    dynamicLhs /= builtinRhs;
    QCOMPARE(TBuiltin(dynamicLhs), builtinLhs);
  }
}

void QDynamicFloatTest::incrementDecrementTestCase()
{
  for(int i = 0; i < iter; ++i)
  {
    //Pre and Post increment
    TBuiltin builtinLhs = dist(e2);
    TDynaFloat dynamicLhs = builtinLhs;
    builtinLhs++;
    dynamicLhs++;
    QCOMPARE(TBuiltin(dynamicLhs), builtinLhs);
    builtinLhs = dist(e2);
    dynamicLhs = builtinLhs;
    ++builtinLhs;
    ++dynamicLhs;
    QCOMPARE(TBuiltin(dynamicLhs), builtinLhs);

    //Pre and Post decrement
    builtinLhs = dist(e2);
    dynamicLhs = builtinLhs;
    builtinLhs--;
    dynamicLhs--;
    QCOMPARE(TBuiltin(dynamicLhs), builtinLhs);
    builtinLhs = dist(e2);
    dynamicLhs = builtinLhs;
    --builtinLhs;
    --dynamicLhs;
    QCOMPARE(TBuiltin(dynamicLhs), builtinLhs);

    //Unary minus
    builtinLhs = dist(e2);
    dynamicLhs = builtinLhs;
    QCOMPARE(TBuiltin(-dynamicLhs), -builtinLhs);
  }
}

QTEST_APPLESS_MAIN(QDynamicFloatTest)

#include "main.moc"
