//
// Created by sba on 19.07.21.
//

#ifndef DUAL_ALGEBRA_KINEMATICS_DUAL_SKEW_H
#define DUAL_ALGEBRA_KINEMATICS_DUAL_SKEW_H

#include "embedded_types/dual_embedded_matrix.h"
#include "screws/unit_line.h"
#include "base/matrix3.h"

class DualSkewProduct;

class DualSkew: public DualEmbeddedMatrix {
public:

    DualSkew(const SkewMatrix &real, const SkewMatrix &dual) noexcept;

    explicit DualSkew(const UnitLine &line) noexcept;

    friend DualSkewProduct;
};

#endif //DUAL_ALGEBRA_KINEMATICS_DUAL_SKEW_H
