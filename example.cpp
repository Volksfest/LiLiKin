//
// Created by sba on 06.07.21.
//

#include <iostream>

#include "dual_number.h"
#include "vector.h"
#include "unit_line.h"
#include "dual_frame.h"

#include "ccc.h"

#include "random.h"

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

ReturnType select(char * input) {
    if ( strncmp("screwed", input, 9) == 0) {
        return screwed();
    } else if ( strncmp("simple_parallel", input, 16) == 0) {
        return simple_parallel();
    } else if ( strncmp("parallel", input, 9) == 0) {
        return parallel();
    } else if ( strncmp("coincide", input, 9) == 0) {
        return coincide();
    } else {
        throw std::logic_error("");
    }
}

void test_scenario(char * input) {
    auto [ccc, goal] = select(input);

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

void test_random() {
    auto a = Random::SampleLine();
    std::cout << "a:" << std::endl << a << std::endl;

    auto b = Random::SampleRelatedLine(a, LineRelation::ANTI_COINCIDE);
    std::cout << "b:" << std::endl << b << std::endl;

    auto [aa,bb,cc] = Random::SampleLineTriplet(LineRelation::PARALLEL, LineRelation::SKEW);
    std::cout << "Triplet a:" << std::endl << aa << "Triplet b:" << std::endl << bb << "Tiplet Ref:" << std::endl << cc << std::endl;
}

int main(int argc, char **argv) {
    try {
        if (argc == 1) {
            throw std::logic_error("");
        }

        if (strncmp("-r", argv[1], 3) == 0) {
            test_random();
        } else {
            test_scenario(argv[1]);
        }
    }
    catch(std::logic_error &e) {
        std::cerr << "Usage: " << std::endl <<
                  argv[0] << " " << "-r" << std::endl <<
                  "or" << std::endl <<
                  argv[0] << " " << "screwed | simple_parallel | parallel | coincide" << std::endl;
    }
}


