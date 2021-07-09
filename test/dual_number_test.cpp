//
// Created by sba on 06.07.21.
//

#include <iostream>
#include "dual_number.h"

#include <gtest/gtest.h>

using namespace DualNumberAlgebra;
using namespace DualNumberAlgebra::literals;

TEST(DualNumerAlgebra, Creation) {
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

TEST(DualNumerAlgebra, Basic) {
    DualNumber a = 1 + 5.5_s;
    DualNumber _a(1,5.5);
    DualNumber b = 3.2 + 8_s;
    DualNumber _b(3.2,8);

    double c = 2.4;

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
}
