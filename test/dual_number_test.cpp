//
// Created by sba on 06.07.21.
//

#include <iostream>
#include <random>
#include <chrono>
#include "dual_number.h"

#include <gtest/gtest.h>

using namespace DualNumberAlgebra;
using namespace DualNumberAlgebra::literals;

#define EXPECT_NEAR_DN(x, y, e) \
EXPECT_NEAR((x).real(), (y).real(), e);\
EXPECT_NEAR((x).dual(), (y).dual(), e)

TEST(DualNumberAlgebra, Creation) { // NOLINT
    DualNumber a_1 = 5;
    DualNumber a_2 = DualNumber(5);
    DualNumber a_3 = DualNumber(5,0);
    DualNumber a_4 = a_2;
    DualNumber a_5(5);
    DualNumber a_6(5,0);
    EXPECT_EQ(a_1, a_2);
    EXPECT_EQ(a_1, a_3);
    EXPECT_EQ(a_1, a_4);
    EXPECT_EQ(a_1, a_5);
    EXPECT_EQ(a_1, a_6);

    DualNumber b_1(0,3);
    DualNumber b_2 = DualNumber(0,3);
    DualNumber b_3 = 3_s;

    EXPECT_EQ(b_1, b_2);
    EXPECT_EQ(b_1, b_3);

    DualNumber c = 5 + 3_s;
    DualNumber mc = -c;
    EXPECT_EQ(mc, -5 - 3_s);
}

TEST(DualNumberAlgebra, Basic) { // NOLINT
    DualNumber a = 1 + 5.5_s;
    DualNumber _a(1,5.5);
    DualNumber b = 3.2 + 8_s;
    DualNumber _b(3.2,8);

    EXPECT_EQ(a, _a);
    EXPECT_EQ(b, _b);

    EXPECT_EQ(a+b, 4.2 + 13.5_s);
    EXPECT_EQ(a-b, -2.2 - 2.5_s);
    EXPECT_EQ(a*b, 3.2 + 25.6_s);

    auto quot = b/a;
    EXPECT_NEAR(quot.real(), 3.2 , 0.01);
    EXPECT_NEAR(quot.dual(), -9.6 , 0.01);

    EXPECT_EQ(-a, -1 - 5.5_s);
    EXPECT_EQ(a.conjugate(), 1 - 5.5_s);
    EXPECT_EQ((-a).dual(), a.conjugate().dual());

    double c = 2.4;
    DualNumber d = c + a;
    EXPECT_EQ(d, 3.4 + 5.5_s);

    DualNumber e = 2 + 2_s;
    EXPECT_EQ(c / e, 1.2 - 1.2_s);
}

TEST(DualNumberAlgebra, ComputeAssign) { //NOLINT
    DualNumber a = 1 + 2_s;
    DualNumber b = 3 - 1_s;
    DualNumber c = 4 + 1_s;

    a += b;
    EXPECT_EQ(a, 4 + 1_s);

    a *= b;
    EXPECT_EQ(a, 12 - 1_s);

    a /= c;
    EXPECT_EQ(a, 3 + -1_s);

    a -=c;
    EXPECT_EQ(a, -1 - 2_s);
}

TEST(DualNumberAlgebra, Advanced) { //NOLINT
    DualNumber a = M_PI_2;
    DualNumber b = DualNumber(0, M_PI_2);

    DualNumber res = sin(a);
    EXPECT_NEAR_DN(res, 1+0_s, 0.01);

    res = cos(a);
    EXPECT_NEAR_DN(res, 0+0_s, 0.01);

    res = sin(b);
    EXPECT_NEAR_DN(res, DualNumber(0,M_PI_2), 0.01);

    res = cos(b);
    EXPECT_NEAR_DN(res, 1+0_s, 0.01);

    res = sin( (a+b) / 2);
    EXPECT_NEAR_DN(res, DualNumber(M_SQRT1_2, M_PI_4 * M_SQRT1_2), 0.01);

    res = cos( (a+b) / 2);
    EXPECT_NEAR_DN(res, DualNumber(M_SQRT1_2, -M_PI_4 * M_SQRT1_2), 0.01);

    res = sqrt( (a+b) / 2);
    EXPECT_NEAR_DN(res, DualNumber(1/M_2_SQRTPI, M_PI_4 * M_2_SQRTPI / 2), 0.01);
}

TEST(DualNumberAlgebra, Norm) { //NOLINT
    DualNumber a = 5 + 3_s;
    DualNumber ma = -a;
    DualNumber b = 5 + 8_s;
    DualNumber c = -5 + 8_s;
    DualNumber d = -2 - 3_s;
    DualNumber e = 4 - 2_s;

    EXPECT_EQ(a.norm(), 5);
    EXPECT_EQ(ma.norm(), 5);
    EXPECT_EQ(b.norm(), 5);
    EXPECT_EQ(c.norm(), 5);

    EXPECT_EQ(d.norm(), 2);
    EXPECT_EQ(e.norm(), 4);

    EXPECT_EQ(a.norm_square(), 25);
    EXPECT_EQ(ma.norm_square(), 25);
    EXPECT_EQ(d.norm_square(), 4);
    EXPECT_EQ(e.norm_square(), 16);
}

TEST(DualNumberAlgebra, Zero) { //NOLINT
    DualNumber z = 0 + 0_s;
    DualNumber zd = 0 + 3_s;
    DualNumber pd = 5 + 0_s;

    DualNumber d = 5 + 3_s;

    EXPECT_TRUE(z.is_zero());
    EXPECT_FALSE(zd.is_zero());
    EXPECT_FALSE(pd.is_zero());

    DualNumber sum = zd + pd - d;
    EXPECT_TRUE(sum.is_zero());
}

TEST(DualNumberAlgebra, Exception) { //NOLINT
    DualNumber z(0, 0);
    DualNumber zd(0, 4);
    DualNumber pd(2, 4);

    EXPECT_NEAR_DN(pd.inverse(), 0.5 - 1_s, 0.01);
    EXPECT_THROW(z.inverse(), std::logic_error);
    EXPECT_THROW(zd.inverse(), std::logic_error);

    EXPECT_NEAR_DN(2 / pd, 1 - 2_s, 0.1);
    EXPECT_THROW(2 / z, std::logic_error);
    EXPECT_THROW(2 / zd, std::logic_error);

    EXPECT_THROW( pd /= zd, std::logic_error);
}

TEST(DualNumberAlgebra, TrigEq) { // NOLINT

}

TEST(DualNumberAlgebra, Random) { // NOLINT
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed); // NOLINT
    std::uniform_real_distribution<double> distribution_to_pi_2(0,M_PI_2);
    std::uniform_real_distribution<double> distribution_circle(-M_PI,M_PI);
    std::uniform_real_distribution<double> distribution_huge(-100,100);
    auto quarter_circle = std::bind(distribution_to_pi_2, generator);
    auto circle = std::bind(distribution_circle, generator);
    auto band = std::bind(distribution_huge, generator);

    DualNumber pos_real(abs(band()), band());
    DualNumber any(quarter_circle(), band());
    EXPECT_NEAR_DN(pos_real, sqrt(pos_real * pos_real), 0.01);
    EXPECT_NEAR_DN(any, acos(cos(any)), 0.01);
    EXPECT_NEAR_DN(any, asin(sin(any)), 0.01);
    EXPECT_NEAR_DN(any, atan(tan(any)), 0.01);

    any = DualNumber(quarter_circle(), band()) - M_PI_2;
    EXPECT_NEAR_DN(pos_real, sqrt(pos_real * pos_real), 0.01);
    EXPECT_NEAR_DN(any, asin(sin(any)), 0.01);
    EXPECT_NEAR_DN(any, atan(tan(any)), 0.01);

    any = DualNumber(quarter_circle(), band()) + M_PI_2;
    EXPECT_NEAR_DN(any, acos(cos(any)), 0.01);

    any = DualNumber(circle(), band());
    EXPECT_NEAR_DN(any, atan2(sin(any), cos(any)), 0.01);
}
