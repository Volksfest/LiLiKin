//
// Created by sba on 20.07.21.
//

#include <iostream>
#include <random>
#include <chrono>

#include "vector.h"
#include "dual_number.h"
#include "screw.h"
#include "unit_line.h"

#include "dual_embedded_matrix.h"
#include "dual_frame.h"
#include "dual_skew.h"
#include "dual_skew_product.h"

#include <gtest/gtest.h>

#define EXPECT_NEAR_DN(x, y, e) \
EXPECT_NEAR((x).real(), (y).real(), e);\
EXPECT_NEAR((x).dual(), (y).dual(), e)

using namespace DualNumberAlgebra;
using namespace DualNumberAlgebra::literals;

TEST(Screws, Creation) { // NOLINT
    Screw a_s(
            DirectionVector(0,0,1),
            MomentVector(0,0,0)
            );

    UnitLine a_l(
            UnitDirectionVector(0,0,1),
            PointVector(0,0,0)
            );

    // Even though a screw is considered to have a transformation, if it doesn't it is practical a line
    EXPECT_EQ(a_s,a_l);

    UnitLine a(
            PointVector(0,0,-1),
            PointVector(0,0,42)
            );

    UnitLine b(
            PointVector(0,1,1),
            PointVector(0,1,2)
            );

    UnitLine c(
            DirectionVector(0,0,1),
            PointVector(0,1,0)
            );

    UnitLine d(
            UnitDirectionVector(0,0,1),
            MomentVector(1,0,0)
            );

    EXPECT_NE(a,b);
    EXPECT_EQ(b,c);
    EXPECT_EQ(b,d);
}

TEST(Screws, Distance) { //NOLINT
    UnitLine z_p(
            UnitDirectionVector(0,0,1),
            PointVector(1,3,0)
            );

    UnitLine z_n(
            UnitDirectionVector(0,0,-1),
            PointVector(3,3,0)
    );

    UnitLine x(
            UnitDirectionVector(1,0,0),
            PointVector(0,2,0)
    );

    UnitLine l(
            DirectionVector(1,1,1),
            PointVector(2,2,0)
            );

    double ang = acos(1/sqrt(3));

    EXPECT_NEAR_DN(z_p.get_distance(z_n), DualNumber(M_PI, -2), 0.000001);
    EXPECT_NEAR_DN(z_p.get_distance(x), DualNumber(M_PI_2, -1), 0.000001);
    EXPECT_NEAR_DN(x.get_distance(l), DualNumber(ang, 0),0.000001);
    EXPECT_NEAR_DN(z_p.get_distance(l), DualNumber(ang, -M_SQRT2), 0.000001);
    EXPECT_NEAR_DN(z_n.get_distance(l), DualNumber(M_PI-ang, 0), 0.000001);
}

TEST(Screws, Orthogonal_Point) { //NOLINT
    UnitLine z_p(
            UnitDirectionVector(0,0,1),
            PointVector(1,3,0)
    );

    UnitLine z_n(
            UnitDirectionVector(0,0,-1),
            PointVector(3,3,0)
    );

    UnitLine x(
            UnitDirectionVector(1,0,0),
            PointVector(0,2,0)
    );

    UnitLine l(
            DirectionVector(1,1,1),
            PointVector(2,2,0)
    );


    auto o_p = z_p.orthogonal(z_n).to_line();
    auto o_d = z_p.orthogonal(l).to_line();
    auto o_x = z_p.orthogonal(x).to_line();

    EXPECT_NEAR_DN(o_p.get_distance(x), DualNumber(M_PI, -1), 0.000001);

    EXPECT_EQ(o_p.n(), DirectionVector(-1,0,0));
    EXPECT_EQ(o_d.n(), DirectionVector(-1,1,0).normal());
    EXPECT_EQ(o_x.n(), DirectionVector(0,1,0).normal());
}

TEST(Screws, Projection) {
    UnitLine l(
            DirectionVector(1,1,0).normal(),
            PointVector(0,0,0)
            );

    PointVector a(0,0,0);
    PointVector b(1,1,1);
    PointVector c(1,0,0);
    PointVector d(1,0,0);

    EXPECT_EQ(l.point_project(a), a);
    EXPECT_EQ(l.point_project(b), PointVector(1,1,0));
    EXPECT_EQ(l.point_project(c), PointVector(0.5,0.5,0));
    EXPECT_EQ(l.orthogonal_plane_projection(PointVector(0,0,0), d), PointVector(0.5, -0.5, 0));
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

    EXPECT_EQ(frame, frame2);
    EXPECT_EQ(skew_invert.skew().screw(), a);
    EXPECT_NEAR_DN(skew_invert.angle(), transform, 0.00001);
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
    EXPECT_NEAR_DN(res, DualNumber(M_PI, -sqrt(2)), 0.00001);

    UnitLine a_p_l(
            UnitDirectionVector(0,0,1),
            PointVector(-1,0,0)
            );
    UnitLine b_p_l(
            UnitDirectionVector(0,0,1),
            PointVector(2,0,0)
            );

    auto res2 = line.acos3(a_p_l, b_p_l);
    EXPECT_NEAR_DN(res2, -M_PI + 0_s, 0.00001);

    auto res3 = line.acos3(line, a_p_l);
    EXPECT_NEAR_DN(res3, 0 + 0_s, 0.000001);

}

TEST(Screws, Intersection) {
    UnitLine line_a(
            UnitDirectionVector(0,1,0),
            PointVector(-1,0,0)
    );

    UnitLine line_b(
            DirectionVector(1,1,0).normal(),
            PointVector(0,0,0)
    );

    UnitLine line_bi(
            DirectionVector(1,1,1).normal(),
            PointVector(0,0,0)
    );

    PointVector intersection = line_a.intersect(line_b);

    EXPECT_EQ(intersection.get(), PointVector(-1,-1,0).get());

    EXPECT_THROW(line_a.intersect(line_bi), std::domain_error);
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


    EXPECT_EQ(a_by_line, a_by_screw);
    EXPECT_EQ(b_by_line, b_by_screw);
    EXPECT_EQ(a_by_line * b_by_line, a_by_screw * b_by_screw);

    //auto sum = skew_a_by_line + skew_b_by_line;
    //auto frame_by_sum = DualFrame(sum);
    //auto skew_by_deconstruction = frame_by_line.constructive_line();
    //ASSERT_EQ(sum, skew_by_deconstruction);
    //ASSERT_EQ(frame_by_sum, a_by_line * b_by_line);
}
