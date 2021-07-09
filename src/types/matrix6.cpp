//
// Created by sba on 07.07.21.
//

#include "types/matrix6.h"
#include "types/matrix3.h"
#include "types/pluecker.h"
#include "types/homogenous.h"

Matrix6::Matrix6(Mat6 data) noexcept {
    this->data = data;
}

Matrix6::Matrix6(const Matrix3 &real, const Matrix3 &dual) noexcept {
    this->data.topLeftCorner(3, 3)     = real.data;
    this->data.topRightCorner(3, 3)     = Eigen::Matrix<double,3,3>::Zero(3, 3);
    this->data.bottomLeftCorner(3, 3)     = dual.data;
    this->data.bottomRightCorner(3, 3)     = real.data;
}

Matrix6::Matrix6(const DualNumberAlgebra::DualNumber &dn) noexcept
    : Matrix6(
            Matrix3(Eigen::Matrix<double,3,3>::Identity(3,3) * dn.real()),
            Matrix3(Eigen::Matrix<double,3,3>::Identity(3,3) * dn.dual())) {}

Matrix6::Matrix6(const Pluecker &pl) noexcept
    : Matrix6(
            SkewMatrix(pl.n()),
            SkewMatrix(pl.m())) {}

Matrix6
Matrix6::operator-() const noexcept {
    return Matrix6(-this->data);
}

Matrix6
Matrix6::operator+(const Matrix6 &rhs) const noexcept {
    return Matrix6(this->data + rhs.data);
}

Matrix6
Matrix6::operator-(const Matrix6 &rhs) const noexcept {
    return Matrix6(this->data - rhs.data);
}

Matrix6
Matrix6::operator*(const Matrix6 &rhs) const noexcept {
    return Matrix6(this->data * rhs.data);
}

Pluecker
Matrix6::operator*(const Pluecker &rhs) const noexcept {
    return Pluecker(this->data * rhs.data);
}

AdjungateMatrix::AdjungateMatrix(const Mat6 &mat) noexcept : Matrix6(mat) {}

AdjungateMatrix::AdjungateMatrix(const Matrix6 &mat) noexcept : Matrix6(mat) {}

AdjungateMatrix::AdjungateMatrix(const RotationMatrix &real, const Matrix3 &dual) noexcept : Matrix6(real,dual) {}

AdjungateMatrix::AdjungateMatrix(const HomogenousMatrix &hom) noexcept
    : Matrix6(
            hom.R(),
            SkewMatrix(hom.p()) * hom.R()
            ) {}

AdjungateMatrix
AdjungateMatrix::operator*(const AdjungateMatrix &rhs) const noexcept {
    return AdjungateMatrix(this->data * rhs.data);
}

Pluecker
AdjungateMatrix::operator*(const Pluecker &rhs) const noexcept {
    return Pluecker(this->data * rhs.data);
}

AdjungateMatrix
AdjungateMatrix::inverse() const noexcept {
    return AdjungateMatrix(this->data.inverse()); // Todo probably faster through transposing
}

RotationMatrix
AdjungateMatrix::R() const noexcept {
    return RotationMatrix(this->data.topLeftCorner(3,3));
}

Matrix3
AdjungateMatrix::pxR() const noexcept {
    return Matrix3(this->data.bottomLeftCorner(3,3));
}

Vector
AdjungateMatrix::p() const noexcept {
    return Vector(SkewMatrix(this->pxR() * this->R().inverse())); // multiple type conversions are needed for safe type conversion
}
