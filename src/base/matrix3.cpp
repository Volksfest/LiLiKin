//
// Created by sba on 07.07.21.
//

#include "base/matrix3.h"
#include "base/vector.h"

#include "util/precision.h"

Matrix3::Matrix3(const Mat3  &data) noexcept: data(data) {}

Matrix3
Matrix3::operator*(double rhs) const noexcept {
    return Matrix3(this->data * rhs);
}

Matrix3
operator*(double lhs, const Matrix3 &rhs) noexcept {
    return rhs * lhs;
}

Matrix3
Matrix3::operator*(const Matrix3 &rhs) const noexcept {
    return Matrix3(this->data * rhs.data);
}

Vector
Matrix3::operator*(const Vector &rhs) const noexcept {
    return Vector(this->data * rhs.get());
}

const Matrix3::Mat3 &
Matrix3::get() const noexcept {
    return data;
}

RotationMatrix::RotationMatrix(double z, double y, double x) noexcept:
    Matrix3(
            (Eigen::AngleAxisd(z, Eigen::Vector3d::UnitZ())*
            Eigen::AngleAxisd(y, Eigen::Vector3d::UnitY())*
            Eigen::AngleAxisd(x, Eigen::Vector3d::UnitX())).matrix()
            ) {}

RotationMatrix::RotationMatrix(const Mat3 &data) noexcept: Matrix3(data) {}

RotationMatrix
RotationMatrix::operator*(const RotationMatrix &rhs) const noexcept {
    return RotationMatrix(this->data * rhs.data);
}

RotationMatrix
RotationMatrix::inverse() const noexcept {
    return RotationMatrix(this->data.transpose());
}

SkewMatrix::SkewMatrix(const Vector & _vector) noexcept: Matrix3(Mat3::Zero(3,3)) {
    auto vector = _vector.get();
    this->data << 0, -vector(2, 0), vector(1, 0),
            vector(2, 0), 0, -vector(0, 0),
            -vector(1, 0), vector(0, 0), 0;
}

SkewMatrix::SkewMatrix(const Matrix3 & rhs) : Matrix3(rhs) {
    auto check = (this->data + this->data.transpose()).array().abs();
    if ( ! (check< Compare::instance().get_precision()).all() ) {
        throw std::domain_error("Matrix is not a skew matrix");
    }
}

RotationMatrix RotationMatrix::RotationFromEigen(const Matrix3::Mat3 &data) {
    auto check = (data - data.transpose()).array().abs();
    if ( ! (check< Compare::instance().get_precision()).all() ) {
        throw std::domain_error("Matrix is not orthogonal");
    }

    if ( !Compare::is_equal(1.0, data.determinant())) {
        throw std::domain_error("Matrix is not orthogonal");
    }

    return RotationMatrix(data);
}
