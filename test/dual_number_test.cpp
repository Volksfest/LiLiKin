//
// Created by sba on 06.07.21.
//

#include <iostream>
#include <random>
#include <chrono>
#include "base/dual_number.h"

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
    EXPECT_NEAR_DN(res, DualNumber(1/M_2_SQRTPIl, M_PI_4 * M_2_SQRTPIl / 2), 0.01);
}

TEST(DualNumberAlgebra, Exception) { //NOLINT
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed); // NOLINT
    std::uniform_real_distribution<double> distribution_huge(-100,100);
    auto band = std::bind(distribution_huge, generator);

    auto zero = DualNumber(0, band());

    auto any = DualNumber(band() + 101, band()); // cannot be zero

    EXPECT_NEAR_DN(any, any.inverse().inverse(), 0.01);
    EXPECT_NEAR_DN(zero, any * (zero / any), 0.01);
    try {
        auto not_existent = zero.inverse();
        FAIL() << not_existent << " is not the inverse of " << zero << "!" << std::endl;
    } catch (std::logic_error &err) {
        EXPECT_EQ(std::string(err.what()), "Cannot invert (a dual number with real part equals to) zero");
    }

    try {
        auto not_existent = any / zero;
        FAIL() << not_existent << " after division should not be existing!" << std::endl;
    } catch (std::logic_error &err) {
        EXPECT_EQ(std::string(err.what()), "Cannot divide by (a real part equal to) zero");
    }
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
