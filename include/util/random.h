//
// Created by sba on 22.11.21.
//

#ifndef DUAL_ALGEBRA_KINEMATICS_RANDOM_H
#define DUAL_ALGEBRA_KINEMATICS_RANDOM_H

#include "base/vector.h"
#include "screws/unit_line.h"
#include "base/matrix3.h"
#include "base/dual_number.h"
#include "embedded_types/dual_frame.h"

namespace Random {
    Vector SampleVector();

    UnitDirectionVector SampleDirection();

    RotationMatrix SampleOrientation();

    UnitLine SampleLine();

    DualNumberAlgebra::DualNumber SampleDualNumber();

    DualFrame SampleFrame();
}

#endif //DUAL_ALGEBRA_KINEMATICS_RANDOM_H
