//
// Created by sba on 08.03.23.
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
EXPECT_TRUE(eq_angle(x,y,e,e))

using namespace DualNumberAlgebra;
using namespace DualNumberAlgebra::literals;

TEST(Acos3, ScrewedScrewed) { // NOLINT
    UnitLine n(
            DirectionVector(0, 0, 1).normal(),
            PointVector(0, 0, 0)
    );

    UnitLine a(
            DirectionVector(0, 0.5, 1.5).normal(),
            PointVector(1, 0, 0)
    );

    DualNumber phi = M_PI_4 + 5_s;

    UnitLine b = DualFrame(DualSkewProduct(n, phi)) * a;
    DualNumber dna = n.get_distance(a);
    DualNumber dnb = n.get_distance(b);

    EXPECT_NEAR_DN(dna, dnb, 1e-8);

    DualNumber c_phi = n.acos3(a,b);
    DualNumber c_inv_phi = n.acos3(b,a);
    EXPECT_NEAR_DN(c_phi, phi, 1e-8);
    EXPECT_NEAR_DN(c_inv_phi, -phi, 1e-8);

    EXPECT_EQ(n.get_relation_to(a), LineRelation::SKEW);
    EXPECT_EQ(a.get_relation_to(b), LineRelation::SKEW);
    }

TEST(Acos3, ScrewedParallel) { // NOLINT
    UnitLine n(
            DirectionVector(0, 1, 1).normal(),
            PointVector(0, 0, 0)
    );

    UnitLine a(
            DirectionVector(0.1, 0.5, 1.5).normal(),
            PointVector(1, 0, 0)
    );

    DualNumber phi = 5_s;

    UnitLine b = DualFrame(DualSkewProduct(n, phi)) * a;
    DualNumber dna = n.get_distance(a);
    DualNumber dnb = n.get_distance(b);

    EXPECT_NEAR_DN(dna, dnb, 1e-8);

    DualNumber c_phi = n.acos3(a,b);
    DualNumber c_inv_phi = n.acos3(b,a);
    EXPECT_NEAR_DN(c_phi, phi, 1e-8);
    EXPECT_NEAR_DN(c_inv_phi, -phi, 1e-8);

    EXPECT_EQ(n.get_relation_to(a), LineRelation::SKEW);
    EXPECT_EQ(a.get_relation_to(b), LineRelation::PARALLEL);
}

TEST(Acos3, ScrewedPlaneParallel) { // NOLINT
    UnitLine n(
            DirectionVector(0, 0, 1).normal(),
            PointVector(0, 0, 0)
    );

    UnitLine a(
            DirectionVector(0, 1, 0).normal(),
            PointVector(1, 0, 0)
    );

    DualNumber phi = M_PI + 3_s;

    UnitLine b = DualFrame(DualSkewProduct(n, phi)) * a;
    DualNumber dna = n.get_distance(a);
    DualNumber dnb = n.get_distance(b);

    EXPECT_NEAR_DN(dna, dnb, 1e-8);

    DualNumber c_phi = n.acos3(a,b);
    DualNumber c_inv_phi = n.acos3(b,a);
    EXPECT_NEAR_DN(c_phi, phi, 1e-8);
    EXPECT_NEAR_DN(c_inv_phi, -phi, 1e-8);

    EXPECT_EQ(n.get_relation_to(a), LineRelation::SKEW);
    EXPECT_EQ(a.get_relation_to(b), LineRelation::ANTI_PARALLEL);
}

TEST(Acos3, ScrewedCoincide) { // NOLINT
    UnitLine n(
            DirectionVector(0, 0, 1).normal(),
            PointVector(0, 0, 0)
    );

    UnitLine a(
            DirectionVector(0, 1, 0).normal(),
            PointVector(1, 0, 0)
    );

    DualNumber phi = 0.0 + 0_s;

    UnitLine b = DualFrame(DualSkewProduct(n, phi)) * a;
    DualNumber dna = n.get_distance(a);
    DualNumber dnb = n.get_distance(b);

    EXPECT_NEAR_DN(dna, dnb, 1e-8);

    DualNumber c_phi = n.acos3(a,b);
    DualNumber c_inv_phi = n.acos3(b,a);
    EXPECT_NEAR_DN(c_phi, phi, 1e-8);
    EXPECT_NEAR_DN(c_inv_phi, -phi, 1e-8);

    EXPECT_EQ(n.get_relation_to(a), LineRelation::SKEW);
    EXPECT_EQ(a.get_relation_to(b), LineRelation::COINCIDE);
}

TEST(Acos3, ParallelParallel) { // NOLINT
    UnitLine n(
            DirectionVector(0, 0, 1).normal(),
            PointVector(0, 0, 0)
    );

    UnitLine a(
            DirectionVector(0, 0, 1).normal(),
            PointVector(1, 0, 0)
    );

    DualNumber phi = 3.0 + 5_s;

    UnitLine b = DualFrame(DualSkewProduct(n, phi)) * a;
    DualNumber dna = n.get_distance(a);
    DualNumber dnb = n.get_distance(b);

    EXPECT_NEAR_DN(dna, dnb, 1e-8);

    DualNumber c_phi = n.acos3(a,b);
    DualNumber c_inv_phi = n.acos3(b,a);
    EXPECT_NEAR(c_phi.real(), phi.real(),1e-8);
    EXPECT_NEAR(c_inv_phi.real(), -phi.real(),1e-8);
    EXPECT_TRUE(std::abs(c_phi.dual() - phi.dual()) > 1e-8);

    EXPECT_EQ(n.get_relation_to(a), LineRelation::PARALLEL);
    EXPECT_EQ(a.get_relation_to(b), LineRelation::PARALLEL);
}

TEST(Acos3, ParallelCoincide) { // NOLINT
    UnitLine n(
            DirectionVector(0, 0, 1).normal(),
            PointVector(0, 0, 0)
    );

    UnitLine a(
            DirectionVector(0, 0, 1).normal(),
            PointVector(1, 0, 0)
    );

    DualNumber phi = 0.0 + 5_s;

    UnitLine b = DualFrame(DualSkewProduct(n, phi)) * a;
    DualNumber dna = n.get_distance(a);
    DualNumber dnb = n.get_distance(b);

    EXPECT_NEAR_DN(dna, dnb, 1e-8);

    DualNumber c_phi = n.acos3(a,b);
    DualNumber c_inv_phi = n.acos3(b,a);
    EXPECT_NEAR(c_phi.real(), phi.real(),1e-8);
    EXPECT_NEAR(c_inv_phi.real(), -phi.real(),1e-8);
    EXPECT_TRUE(std::abs(c_phi.dual() - phi.dual()) > 1e-8);

    EXPECT_EQ(n.get_relation_to(a), LineRelation::PARALLEL);
    EXPECT_EQ(a.get_relation_to(b), LineRelation::COINCIDE);
}

TEST(Acos3, CoincideCoincide) { // NOLINT
    UnitLine n(
            DirectionVector(0, 0, 1).normal(),
            PointVector(0, 0, 0)
    );

    UnitLine a = n;

    DualNumber phi = 1.0 + 5_s;

    UnitLine b = DualFrame(DualSkewProduct(n, phi)) * a;
    DualNumber dna = n.get_distance(a);
    DualNumber dnb = n.get_distance(b);

    EXPECT_NEAR_DN(dna, dnb, 1e-8);

    DualNumber c_phi = n.acos3(a,b);
    EXPECT_NEAR_DN(c_phi, 0 + 0_s, 1e-8);

    EXPECT_EQ(n.get_relation_to(a), LineRelation::COINCIDE);
    EXPECT_EQ(a.get_relation_to(b), LineRelation::COINCIDE);
}