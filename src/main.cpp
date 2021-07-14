//
// Created by sba on 06.07.21.
//

#include <iostream>
#include "dual_number.h"

#include "types.h"
#include "ccc.h"

using namespace DualNumberAlgebra;

DualNumber toDeg(const DualNumber &a) {
    return DualNumber(a.real()/M_PI*180.0, a.dual());
}

int main() {
    Pluecker a(
            DirectionVector(0.0, 0.0, 1.0),
            PointVector(0.0, 0.0, 0.0)
    );

    Pluecker b(
            DirectionVector(1.0, 0.0, 0.0),
            PointVector(0.0, 0.0, 0.0)
    );

    Pluecker c(
            DirectionVector(0.0, 0.0, 1.0),
            PointVector(0.0, 0.0, 0.0)
    );

    HomogenousMatrix zero_posture_hom(
            RotationMatrix(0,0,0),
            Vector(0,0,4)
    );

    AdjungateMatrix zero_posture(zero_posture_hom);

    HomogenousMatrix to_hom(
            RotationMatrix(0,0,0),
            Vector(-5,5,8)
    );

    AdjungateMatrix to(to_hom);

    CCCMechanism ccc(a,b,c, zero_posture);

    auto solutions = ccc.inverse(to);

    std::cout << "To Pose: " << std::endl << to_hom << std::endl;

    for (auto solution : solutions) {
        std::cout << std::endl <<
                  toDeg(solution.phi_1) << std::endl <<
                  toDeg(solution.phi_2) << std::endl <<
                  toDeg(solution.phi_3) << std::endl;

        std::cout << std::endl <<
                  "Yields to: " << std::endl << HomogenousMatrix(ccc.forward(solution)) << std::endl;
    }
}
