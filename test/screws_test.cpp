//
// Created by sba on 20.07.21.
//

#include <iostream>
#include <random>
#include <chrono>

#include "base/vector.h"
#include "screws/screw.h"
#include "screws/line.h"
#include "screws/unit_screw.h"
#include "screws/unit_line.h"

#include <gtest/gtest.h>

#define EXPECT_NEAR_DN(x, y, e) \
EXPECT_NEAR((x).real(), (y).real(), e);\
EXPECT_NEAR((x).dual(), (y).dual(), e)

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