//
// Created by sba on 20.07.21.
//

#include <iostream>
#include <random>
#include <chrono>

#include "base/vector.h"
#include "base/dual_number.h"
#include "screws/screw.h"
#include "screws/line.h"
#include "screws/unit_screw.h"
#include "screws/unit_line.h"

#include "embedded_types/dual_embedded_matrix.h"
#include "embedded_types/dual_frame.h"
#include "embedded_types/dual_skew.h"
#include "embedded_types/dual_skew_product.h"

#include <gtest/gtest.h>

#define EXPECT_NEAR_DN(x, y, e) \
EXPECT_NEAR((x).real(), (y).real(), e);\
EXPECT_NEAR((x).dual(), (y).dual(), e)

using namespace DualNumberAlgebra;
using namespace DualNumberAlgebra::literals;

TEST(Screws, Creation) { // NOLINT
    Screw a(
            DirectionVector(0,0,1),
            MomentVector(0,0,0)
            );
    Line b(
            DirectionVector(0,0,1),
            PointVector(0,0,1)
            );

    UnitScrew c(
            UnitDirectionVector(0,0,1),
            MomentVector(0,0,0)
            );
    UnitLine d(
            PointVector(0,0,-1),
            PointVector(0,0,42)
            );

    ASSERT_EQ(a,b);
    ASSERT_EQ(c,d);
    ASSERT_NE(a,c);
}

TEST(Screws, Transformation) { // NOLINT
    UnitLine unit_line(
            UnitDirectionVector(0,0,1),
            PointVector(0,0,0)
            );

    auto unit_screw = unit_line.translate(1.5);
    auto line = unit_line.rotate(M_PI_4);

    auto screw = unit_line.transform(DualNumber(M_PI_4, 1.5));
    auto line_to_screw = line.translate(1.5);
    auto unit_screw_to_screw = unit_screw.rotate(M_PI_4);

    ASSERT_EQ(screw, line_to_screw);
    ASSERT_EQ(screw, unit_screw_to_screw);
    ASSERT_EQ(line_to_screw, unit_screw_to_screw);
}

TEST(Screws, Distance) { //NOLINT
    UnitLine a(
            UnitDirectionVector(0,0,1),
            PointVector(0,1,0)
            );

    UnitLine b(
            UnitDirectionVector(0,0,1),
            PointVector(1,0,0)
            );

    UnitLine bi(
            UnitDirectionVector(0,0,-1),
            PointVector(1,0,0)
            );

    UnitLine c(
            UnitDirectionVector(1,0,0),
            PointVector(0,5,0)
            );

    EXPECT_NEAR_DN(a.get_distance(b), DualNumber(0, M_SQRT2), 0.000001);
    EXPECT_NEAR_DN(b.get_distance(bi), DualNumber(M_PI, 0), 0.000001);
    EXPECT_NEAR_DN(a.get_distance(c), DualNumber(M_PI_2, 4),0.000001);
}

TEST(Screws, Inverse) { //NOLINT
    UnitLine a(
            UnitDirectionVector(0,0,1),
            PointVector(1,1,0)
            );

    DualNumber transform(M_PI_4, 5);

    auto skew = DualSkewProduct(a, transform);
    auto frame = DualFrame(skew);
    auto skew_invert = frame.constructive_line();
    auto frame2 = DualFrame(skew_invert);
    ASSERT_EQ(frame, frame2);

}

TEST(Screws, Acos3) { //NOLINT
    UnitLine line(
            UnitDirectionVector(0,0,1),
            PointVector(0,0,0)
            );

    UnitLine a(
            UnitDirectionVector(1,0,0),
            PointVector(0,2,1)
            );

    UnitLine b1(
            UnitDirectionVector(0,1,0),
            PointVector(1,0,1)
            );

    UnitLine b2(
            UnitDirectionVector(0,1,0),
            PointVector(1,0,2)
            );

    EXPECT_NEAR_DN(line.acos3(a,b1), M_PI_2 + 0_s, 0.00001);
    EXPECT_NEAR_DN(line.acos3(a,b2), M_PI_2 + 1_s, 0.00001);
    EXPECT_NEAR_DN(line.acos3(b1,b2), 1_s, 0.00001);

    UnitLine line2(
            DirectionVector(1,1,0).normal(),
            PointVector(1,1,1)
            );

    auto res = line2.acos3(a,b1);
    EXPECT_NEAR_DN(res, DualNumber(M_PI, sqrt(2)), 0.00001);

    UnitLine a_p_l(
            UnitDirectionVector(0,0,1),
            PointVector(-1,0,0)
            );
    UnitLine b_p_l(
            UnitDirectionVector(0,0,1),
            PointVector(2,0,0)
            );

    auto res2 = line.acos3(a_p_l, b_p_l);
    EXPECT_NEAR_DN(res2, M_PI + 0_s, 0.00001);

}

TEST(Screws, LieAlgebra) { //NOLINT
    UnitLine line_a(
            UnitDirectionVector(0,0,1),
            PointVector(0,0,0)
            );
    UnitLine line_b(
            UnitDirectionVector(0,1,0),
            PointVector(1,0,1)
            );

    DualNumber value_a(M_PI_4, 2.5);
    DualNumber value_b(M_PI_2, 3);

    auto screw_a = line_a.transform(value_a);
    auto screw_b = line_b.transform(value_b);

    auto skew_a_by_line = DualSkewProduct(line_a, value_a);
    auto skew_b_by_line = DualSkewProduct(line_b, value_b);

    auto a_by_line = DualFrame(skew_a_by_line);
    auto b_by_line = DualFrame(skew_b_by_line);

    auto a_by_screw = DualFrame(DualSkewProduct(screw_a));
    auto b_by_screw = DualFrame(DualSkewProduct(screw_b));

    auto frame_by_line = a_by_line * b_by_line;


    ASSERT_EQ(a_by_line, a_by_screw);
    ASSERT_EQ(b_by_line, b_by_screw);
    ASSERT_EQ(a_by_line * b_by_line, a_by_screw * b_by_screw);

    //auto sum = skew_a_by_line + skew_b_by_line;
    //auto frame_by_sum = DualFrame(sum);
    //auto skew_by_deconstruction = frame_by_line.constructive_line();
    //ASSERT_EQ(sum, skew_by_deconstruction);
    //ASSERT_EQ(frame_by_sum, a_by_line * b_by_line);
}
