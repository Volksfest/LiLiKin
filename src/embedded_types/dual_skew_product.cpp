//
// Created by sba on 20.07.21.
//

#include "embedded_types/dual_embedded_matrix.h"
#include "embedded_types/dual_skew_product.h"
#include "screws/screw.h"
#include "base/dual_number.h"
#include "base/vector.h"

DualSkewProduct::DualSkewProduct(const Screw &screw) noexcept:
    _skew(screw.to_line()), _angle(screw.norm()) {}

DualSkewProduct::DualSkewProduct(const UnitLine &screw, const DualNumberAlgebra::DualNumber &argument) noexcept:
    _skew(screw), _angle(argument) {}

DualSkew DualSkewProduct::skew() const noexcept {
    return this->_skew;
}
DualNumberAlgebra::DualNumber DualSkewProduct::angle() const noexcept {
    return this->_angle;
}

DualSkewProduct DualSkewProduct::operator+(const DualSkewProduct &rhs) const noexcept {
    //Eigen::Matrix<double,6,6> sum = DualEmbeddedMatrix(this->_angle).data * this->_skew.data + DualEmbeddedMatrix(rhs._angle).data * rhs._skew.data;

    Eigen::Matrix<double, 6,6> weighted_lhs_raw = DualEmbeddedMatrix(this->_angle).get() * this->_skew.get();
    Eigen::Matrix<double, 6,6> weighted_rhs_raw = DualEmbeddedMatrix(rhs._angle).get() * rhs._skew.get();

    Eigen::Matrix<double, 6,6> sum = weighted_lhs_raw * weighted_rhs_raw - weighted_rhs_raw * weighted_lhs_raw; //Lie Bracket

    // To create type safety, sum will be decomposed into the skews itself and then be rebuild.
    auto n_skew = SkewMatrix(Matrix3(sum.topLeftCorner(3,3)));
    auto m_skew = SkewMatrix(Matrix3(sum.bottomLeftCorner(3,3)));
    // this step is allthough still necessary as DualSkewProduct contains the unit line and thus normalizes the screw
    // The norm itself is again the dual angle and saved as a second variable
    auto screw = Screw(DirectionVector(Vector(n_skew)), MomentVector(Vector(m_skew)));

    return {screw};
}

bool operator==(const DualSkewProduct &lhs, const DualSkewProduct &rhs) noexcept {
    return lhs._angle == rhs._angle && lhs._skew == rhs._skew;
}
