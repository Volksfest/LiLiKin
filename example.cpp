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
using ReturnType = std::pair<CCCMechanism, DualFrame>;

DualNumber toDeg(const DualNumber &a) {
    return DualNumber(a.real()/M_PI*180.0, a.dual());
}

ReturnType coincide() {
    UnitLine a(
            UnitDirectionVector(0.0, 0.0, 1.0),
            PointVector(0.0, 0.0, 0.0)
    );
    UnitLine b(
            UnitDirectionVector(1.0, 0.0, 0.0),
            PointVector(0.0, -1.0, 0.0)
    );
    UnitLine c(
            DirectionVector(0.0, 1.0, 0.0),
            PointVector(10.0, 0.0, 8.0)
    );
    DualFrame zero_posture(
            RotationMatrix(M_PI_2,M_PI_2,0),
            PointVector(10,5,7)
    );

    CCCMechanism ccc(a,b,c, zero_posture);
    DualFrame goal(
            RotationMatrix(0,0,0),
            PointVector(15,0,0)
    );

    return {ccc, goal};
}

ReturnType parallel() {
    UnitLine a(
            UnitDirectionVector(0.0, 0.0, 1.0),
            PointVector(0.0, 0.0, 0.0)
    );
    UnitLine b(
            UnitDirectionVector(1.0, 0.0, 0.0),
            PointVector(0.0, -1.0, 0.0)
    );
    UnitLine c(
            DirectionVector(0.0, -1.0, 1.0),
            PointVector(1.0, 0.0, 0.0)
    );
    DualFrame zero_posture(
            RotationMatrix(0,0,0),
            PointVector(0,-1,4)
    );
    CCCMechanism ccc(a,b,c, zero_posture);
    DualFrame goal(
            RotationMatrix(0,0,-M_PI_4),
            PointVector(4,4,0)
    );
    return {ccc, goal};
}

ReturnType simple_parallel() {
    UnitLine a(
            UnitDirectionVector(0.0, 0.0, 1.0),
            PointVector(0.0, 0.0, 0.0)
    );
    UnitLine b(
            UnitDirectionVector(1.0, 0.0, 0.0),
            PointVector(0.0, -1.0, 0.0)
    );
    UnitLine c(
            DirectionVector(0.0, 0, 1.0),
            PointVector(2.0, 1.0, 0.0)
    );
    DualFrame zero_posture(
            RotationMatrix(0,0,0),
            PointVector(2,0,0)
    );
    CCCMechanism ccc(a,b,c, zero_posture);
    DualFrame goal(
            RotationMatrix(0,0,0),
            PointVector(4,4,4)
    );
    return {ccc, goal};
}

ReturnType screwed() {
    UnitLine a(
            UnitDirectionVector(0.0, 0.0, 1.0),
            PointVector(0.0, 0.0, 0.0)
    );
    UnitLine b(
            UnitDirectionVector(1.0, 0.0, 0.0),
            PointVector(0.0, -1.0, 0.0)
    );
    UnitLine c(
            DirectionVector(0.0, 1.0, 0.0),
            PointVector(1.0, 0.0, 1.0)
    );
    DualFrame zero_posture(
            RotationMatrix(0,0,0),
            PointVector(0,-1,4)
    );
    CCCMechanism ccc(a,b,c, zero_posture);
    DualFrame goal(
            RotationMatrix(M_PI_2,0,-M_PI_4),
            PointVector(4,4,0)
    );
    return {ccc, goal};
}

int main() {
    auto [ccc, goal] = coincide();

    auto solutions = ccc.inverse(goal);

    std::cout << "To Pose: " << std::endl << goal << std::endl;

    for (auto solution : solutions) {
        std::cout << std::endl <<
                  solution.phi_1 << std::endl <<
                  solution.phi_2 << std::endl <<
                  solution.phi_3 << std::endl;
        std::cout << "yields to:" << std::endl << ccc.forward(solution) << std::endl;
    }
}


