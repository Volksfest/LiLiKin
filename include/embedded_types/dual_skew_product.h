//
// Created by sba on 20.07.21.
//

#ifndef DUAL_ALGEBRA_KINEMATICS_DUAL_SKEW_PRODUCT_H
#define DUAL_ALGEBRA_KINEMATICS_DUAL_SKEW_PRODUCT_H

#include "embedded_types/dual_skew.h"
#include "base/dual_number.h"
#include "screws/screw.h"

class DualSkewProduct {
private:
    DualSkew _skew;
    DualNumberAlgebra::DualNumber _angle;
public:
    DualSkewProduct() = delete;
    DualSkewProduct(const Screw &screw) noexcept;
    DualSkewProduct(const UnitLine &screw, const DualNumberAlgebra::DualNumber &argument) noexcept;

    DualSkew skew() const noexcept;
    DualNumberAlgebra::DualNumber angle() const noexcept;

    // this does not work as it is.
    // it was actually just a try in the first place
    DualSkewProduct operator+(const DualSkewProduct &rhs) const noexcept;

    friend bool operator==(const DualSkewProduct &lhs, const DualSkewProduct &rhs) noexcept;
};

#endif //DUAL_ALGEBRA_KINEMATICS_DUAL_SKEW_PRODUCT_H
