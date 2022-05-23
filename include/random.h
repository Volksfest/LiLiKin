//
// Created by sba on 22.11.21.
//

#ifndef DUAL_ALGEBRA_KINEMATICS_RANDOM_H
#define DUAL_ALGEBRA_KINEMATICS_RANDOM_H

#include "vector.h"
#include "unit_line.h"
#include "matrix3.h"
#include "dual_number.h"
#include "dual_frame.h"

namespace Random {
    Vector SampleVector();

    UnitDirectionVector SampleDirection();

    RotationMatrix SampleOrientation();

    UnitLine SampleLine();

    DualNumberAlgebra::DualNumber SampleDualNumber();

    DualFrame SampleFrame();

    PointVector SamplePoint(double cube_size = 1);

    UnitLine SampleRelatedLine(const UnitLine &a, LineRelation relation);

    std::pair<UnitLine, UnitLine> SampleLinePair(LineRelation relation);

    std::tuple<UnitLine, UnitLine, UnitLine> SampleLineTriplet(LineRelation relation_a_to_b, LineRelation relation_a_to_ref);
}

#endif //DUAL_ALGEBRA_KINEMATICS_RANDOM_H
