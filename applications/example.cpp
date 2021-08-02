//
// Created by sba on 06.07.21.
//

#include <iostream>

#include "base/dual_number.h"
#include "base/vector.h"
#include "screws/unit_line.h"
#include "embedded_types/dual_frame.h"

#include "ccc.h"

using namespace DualNumberAlgebra;

DualNumber toDeg(const DualNumber &a) {
    return DualNumber(a.real()/M_PI*180.0, a.dual());
}

int main() {
    UnitLine a(
            UnitDirectionVector(0.0, 0.0, 1.0),
            PointVector(0.0, 0.0, 0.0)
    );

    UnitLine b(
            UnitDirectionVector(1.0, 0.0, 0.0),
            PointVector(0.0, 0.0, 0.0)
    );

    UnitLine c(
            UnitDirectionVector(0.0, 0.0, 1.0),
            PointVector(0.0, 0.0, 0.0)
    );

    DualFrame zero_posture(
            RotationMatrix(0,0,0),
            PointVector(0,0,4)
    );

    DualFrame goal(
            RotationMatrix(0,0,0),
            PointVector(-5,5,8)
    );

    CCCMechanism ccc(a,b,c, zero_posture);

    auto solutions = ccc.inverse(goal);

    std::cout << "To Pose: " << std::endl << goal << std::endl;

    for (auto solution : solutions) {
        std::cout << std::endl <<
                  toDeg(solution.phi_1) << std::endl <<
                  toDeg(solution.phi_2) << std::endl <<
                  toDeg(solution.phi_3) << std::endl;

        std::cout << std::endl <<
                  "Yields to: " << std::endl << ccc.forward(solution) << std::endl;
    }
}


