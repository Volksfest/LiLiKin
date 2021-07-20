//
// Created by sba on 20.07.21.
//

#include "embedded_types/dual_embedded_matrix.h"
#include "embedded_types/dual_skew_product.h"
#include "screws/screw.h"
#include "screws/line.h"
#include "base/dual_number.h"
#include "base/vector.h"

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
    auto sum = DualEmbeddedMatrix(this->_angle).data * this->_skew.data + DualEmbeddedMatrix(rhs._angle).data * rhs._skew.data;

    auto n_skew = SkewMatrix(Matrix3(sum.topLeftCorner(3,3)));
    auto m_skew = SkewMatrix(Matrix3(sum.bottomLeftCorner(3,3)));
    auto screw = Screw(DirectionVector(Vector(n_skew)), MomentVector(Vector(m_skew)));
    auto line = screw.align().normalize();
    auto norm = screw.norm();
    return {line, norm};
}