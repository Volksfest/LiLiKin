//
// Created by sba on 07.07.21.
//

#include "types/matrix3.h"
#include "types/pluecker.h"

Matrix3::Matrix3(Mat3  data) noexcept {
    this->data = data;
}

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
    return Vector(this->data * rhs.data);
}

RotationMatrix::RotationMatrix(double z, double y, double x) noexcept {
    Eigen::AngleAxisd rollAngle(z, Eigen::Vector3d::UnitZ());
    Eigen::AngleAxisd yawAngle(y, Eigen::Vector3d::UnitY());
    Eigen::AngleAxisd pitchAngle(x, Eigen::Vector3d::UnitX());

    Eigen::Quaternion<double> q = rollAngle * yawAngle * pitchAngle;

    this->data = q.matrix();
}

RotationMatrix::RotationMatrix(Mat3 data) noexcept {
    this->data = data;
}

RotationMatrix
RotationMatrix::operator*(const RotationMatrix &rhs) const noexcept {
    return RotationMatrix(this->data * rhs.data);
}

RotationMatrix
RotationMatrix::inverse() const noexcept {
    return RotationMatrix(this->data.transpose());
}

SkewMatrix::SkewMatrix(const Vector & vector) noexcept {
    this->data << 0, -vector.data(2, 0), vector.data(1, 0),
            vector.data(2, 0), 0, -vector.data(0, 0),
            -vector.data(1, 0), vector.data(0, 0), 0;
}

SkewMatrix::SkewMatrix(const Matrix3 & rhs) : Matrix3(rhs) {
    if ( (this->data + this->data.transpose()) != Mat3::Zero(3,3) ) {
        throw std::logic_error("Matrix is not a skew matrix");
    }
}
