//
// Created by sba on 20.07.21.
//

#include "embedded_types/dual_skew_product.h"
#include "screws/screw.h"
#include "screws/line.h"
#include "base/dual_number.h"

DualSkewProduct::DualSkewProduct(const Screw &screw) noexcept:
    _skew(screw.align().normalize()), _angle(screw.norm()) {}

DualSkewProduct::DualSkewProduct(const UnitLine &screw, const DualNumberAlgebra::DualNumber &argument) noexcept:
    _skew(screw), _angle(argument) {}

DualSkew DualSkewProduct::skew() const noexcept {
    return this->_skew;
}
DualNumberAlgebra::DualNumber DualSkewProduct::angle() const noexcept {
    return this->_angle;
}

DualSkewProduct DualSkewProduct::operator+(const DualSkewProduct &rhs) const noexcept {
    return *this; //TODO change!
}