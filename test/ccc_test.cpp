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

using namespace DualNumberAlgebra;
using namespace DualNumberAlgebra::literals;

class Mechanism : public testing::Test {
protected:
    void SetUp() override {
    UnitLine z(
            UnitDirectionVector(0,0,1),
            PointVector(0,0,0)
    );

    UnitLine x(
            UnitDirectionVector(1,0,0),
            PointVector(0,0,0)
    );

    DualFrame su_zp(
            RotationMatrix(0,0,0),
            PointVector(0,0,0)
    );

    SU = std::make_unique<CCCMechanism>(z,x,z,su_zp);






    }

    std::vector<DualFrame> frames;

    std::unique_ptr<CCCMechanism> SU;
    std::unique_ptr<CCCMechanism> orthogonal;
    std::unique_ptr<CCCMechanism> arbitrary;
};

// Helper function

CCCMechanism create_SU() {
    UnitLine z(
            UnitDirectionVector(0,0,1),
            PointVector(0,0,0)
    );

    UnitLine x(
            UnitDirectionVector(1,0,0),
            PointVector(0,0,0)
    );

    DualFrame su_zp(
            RotationMatrix(0,0,0),
            PointVector(0,0,0)
    );

    return CCCMechanism(z,x,z,su_zp);
}

void systematic(CCCMechanism &m) {
    for (int rz_i = -2; rz_i < 2; rz_i++) {
        for (int ry_i = -2; ry_i < 2; ry_i++) {
            for (int rx_i = -2; rx_i < 2; rx_i++) {
                for (int z_i = -2; z_i < 2; z_i++) {
                    for (int y_i = -2; y_i < 2; y_i++) {
                        for (int x_i = -2; x_i < 2; x_i++) {
                            DualFrame frame(
                                    RotationMatrix(rz_i * M_PI_4, ry_i * M_PI_4, rx_i * M_PI_4),
                                    PointVector(x_i, y_i, z_i)
                            );
                            auto configs = m.inverse(frame);

                            for (const auto &config : configs) {
                                auto fk = m.forward(config);
                                EXPECT_EQ(fk, frame) << "(p: " << x_i << ", " << y_i << ", " << z_i << ")" << std::endl
                                                     << "(R: " << rz_i << ", " << ry_i << ", " << rx_i << ")" << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }
}

TEST(Mechanism, SU_Fixed_Frames) { // NOLINT
    CCCMechanism su = create_SU();

    std::vector<DualFrame> frames = {
            DualFrame(RotationMatrix(0,0,0), PointVector(0,0,0)),
            DualFrame(RotationMatrix(0,0,0), PointVector(0,0,-5)),
            DualFrame(RotationMatrix(2,0,0), PointVector(0,0,-5)),
            DualFrame(RotationMatrix(M_PI,0,0), PointVector(0,0,-5)),
            DualFrame(RotationMatrix(M_PI_2,0,0), PointVector(0,0,-5)),
            DualFrame(RotationMatrix(0,M_PI_2,0), PointVector(0,0,-5)),
            DualFrame(RotationMatrix(M_PI_2,M_PI_2,0), PointVector(0,0,-5)),
            DualFrame(RotationMatrix(0,M_PI,0), PointVector(0,0,-5)),
            DualFrame(RotationMatrix(M_PI,M_PI,0), PointVector(0,0,-5)),
            DualFrame(RotationMatrix(1,M_PI,2), PointVector(0,0,-5)),
            DualFrame(RotationMatrix(0,0,0), PointVector(1,0,-5)),
            DualFrame(RotationMatrix(2,0,0), PointVector(0,1,-5)),
            DualFrame(RotationMatrix(M_PI,0,0), PointVector(2,0,-5)),
            DualFrame(RotationMatrix(M_PI_2,0,0), PointVector(0,2,-5)),
            DualFrame(RotationMatrix(0,M_PI_2,0), PointVector(3,0,-5)),
            DualFrame(RotationMatrix(M_PI_2,M_PI_2,0), PointVector(0,3,-5)),
            DualFrame(RotationMatrix(0,M_PI,0), PointVector(4,0,-5)),
            DualFrame(RotationMatrix(M_PI,M_PI,0), PointVector(0,4,-5)),
            DualFrame(RotationMatrix(1,M_PI,2), PointVector(5,5,-5)),
            DualFrame(RotationMatrix(M_PI_4, 0, M_PI_4), PointVector(1, 0, 0))
    };

    for (const auto &frame : frames) {
        auto configs = su.inverse(frame);

        for (const auto &config : configs) {
            auto fk = su.forward(config);
            EXPECT_EQ(fk, frame);
        }
    }
}

TEST(Mechanism, SU_Systematic_Frames) { // NOLINT
    auto mechanism = create_SU();
    systematic(mechanism);
}

TEST(Mechanism, Orthogonal_Systematic_Frames) {
    UnitLine a(
            DirectionVector(1,0,1).normal(),
            PointVector(0,0,0)
    );

    UnitLine b(
            DirectionVector(0,1,0).normal(),
            PointVector(1,0,0)
    );

    UnitLine c(
            DirectionVector(1,0,0).normal(),
            PointVector(0,-4,1)
    );

    DualFrame zp(
            RotationMatrix(1 * M_PI_4, -1 * M_PI_4, 3 * M_PI_4),
            PointVector(-2,0,4)
    );

    CCCMechanism orthogonal(a, b, c, zp);
    systematic(orthogonal);
}

TEST(Mechanism, Orthogonal_Enforced_Frames) { // NOLINT
    UnitLine a(
            DirectionVector(1,0,1).normal(),
            PointVector(0,0,0)
    );

    UnitLine b(
            DirectionVector(0,1,0).normal(),
            PointVector(1,0,0)
    );

    UnitLine c(
            DirectionVector(1,0,0).normal(),
            PointVector(0,2,M_SQRT1_2)
    );

    DualFrame zp(
            RotationMatrix(1 * M_PI_4, -1 * M_PI_4, 3 * M_PI_4),
            PointVector(-2,0,4)
    );

    CCCMechanism orthogonal(a,b,c,zp);

    std::vector<Configuration> configs = {
            {0, -2_s, 0}, // Intersecting
            {0, -M_PI_4, 0}, // Parallel
            {0, -M_PI_4-2_s,0} // Coinciding
    };

    for (auto config : configs) {
        auto frame = orthogonal.forward(config);
        auto inverse_configs = orthogonal.inverse(frame);
        for (const auto &inverse_config : inverse_configs) {
            auto fk = orthogonal.forward(config);
            EXPECT_EQ(fk, frame);
        }
    }
}
