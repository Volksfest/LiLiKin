//
// Created by sba on 19.07.21.
//

#include "dual_embedded_matrix.h"

#include "screw.h"

#include "matrix3.h"
#include "dual_number.h"

DualEmbeddedMatrix::DualEmbeddedMatrix(const Mat6 &data) noexcept: data(data) {}

DualEmbeddedMatrix::DualEmbeddedMatrix(const Matrix3 &real, const Matrix3 &dual) noexcept {
    this->data.topLeftCorner(3, 3)     = real.get();
    this->data.topRightCorner(3, 3)     = Eigen::Matrix<double,3,3>::Zero(3, 3);
    this->data.bottomLeftCorner(3, 3)     = dual.get();
    this->data.bottomRightCorner(3, 3)     = real.get();
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

const DualEmbeddedMatrix::Mat6 &
DualEmbeddedMatrix::get() const noexcept {
    return data;
}
