//
// Created by sba on 27.07.21.
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

#include "ccc.h"

#include <gtest/gtest.h>

#define EXPECT_NEAR_DN(x, y, e) \
EXPECT_NEAR((x).real(), (y).real(), e);\
EXPECT_NEAR((x).dual(), (y).dual(), e)

using namespace DualNumberAlgebra;
using namespace DualNumberAlgebra::literals;

TEST(Mechanism, SU_Fixed_Frames) { // NOLINT

    UnitLine a(
            UnitDirectionVector(0,0,1),
            PointVector(0,0,0)
            );

    UnitLine b(
            UnitDirectionVector(1,0,0),
            PointVector(0,0,0)
    );

    UnitLine c(
            UnitDirectionVector(0,0,1),
            PointVector(0,0,0)
    );

    DualFrame zp(RotationMatrix(0,0,0), PointVector(0,0,0));

    CCCMechanism m(a,b,c,zp);

    std::vector<std::pair<DualFrame,int>> frames = {
            std::make_pair(DualFrame(RotationMatrix(0,0,0), PointVector(0,0,0)), 2)
    };

    for (const auto &[frame,size] : frames) {
        auto configs = m.inverse(frame);
        ASSERT_EQ(configs.size(), size);
        for (const auto &config : configs) {
            auto fk = m.forward(config);
            ASSERT_EQ(fk, frame);
        }
    }
}
