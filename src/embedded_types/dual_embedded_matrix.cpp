//
// Created by sba on 19.07.21.
//

#include "embedded_types/dual_embedded_matrix.h"

#include "screws/screw.h"

#include "base/matrix3.h"
#include "base/dual_number.h"

DualEmbeddedMatrix::DualEmbeddedMatrix(const Mat6 &data) noexcept: data(data) {}

DualEmbeddedMatrix::DualEmbeddedMatrix(const Matrix3 &real, const Matrix3 &dual) noexcept {
    this->data.topLeftCorner(3, 3)     = real.data;
    this->data.topRightCorner(3, 3)     = Eigen::Matrix<double,3,3>::Zero(3, 3);
    this->data.bottomLeftCorner(3, 3)     = dual.data;
    this->data.bottomRightCorner(3, 3)     = real.data;
}

DualEmbeddedMatrix::DualEmbeddedMatrix(const DualNumberAlgebra::DualNumber &dn) noexcept
        : DualEmbeddedMatrix(
        Matrix3(Eigen::Matrix<double,3,3>::Identity(3,3) * dn.real()),
        Matrix3(Eigen::Matrix<double,3,3>::Identity(3,3) * dn.dual())) {}

DualEmbeddedMatrix DualEmbeddedMatrix::operator-() const noexcept {
    return DualEmbeddedMatrix(-this->data);
}

DualEmbeddedMatrix DualEmbeddedMatrix::operator+(const DualEmbeddedMatrix &rhs) const noexcept {
    return DualEmbeddedMatrix(this->data + rhs.data);
}

DualEmbeddedMatrix DualEmbeddedMatrix::operator-(const DualEmbeddedMatrix &rhs) const noexcept {
    return DualEmbeddedMatrix(this->data - rhs.data);
}

DualEmbeddedMatrix DualEmbeddedMatrix::operator*(const DualEmbeddedMatrix &rhs) const noexcept {
    return DualEmbeddedMatrix(this->data * rhs.data);
}
