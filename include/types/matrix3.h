//
// Created by sba on 06.07.21.
//

#ifndef DAK_MATRIX3_H
#define DAK_MATRIX3_H

#include <eigen3/Eigen/Eigen>

class Vector;
class RotationMatrix;
class SkewMatrix;

class Matrix6;
class AdjungateMatrix;
class HomogenousMatrix;

class Matrix3 {
private:
    using Mat3 = Eigen::Matrix<double,3,3>;

    Eigen::Matrix<double,3,3> data;

    Matrix3() = default;
    explicit Matrix3(Mat3  data) noexcept;

public:
    Matrix3 operator*(double rhs) const noexcept;
    Matrix3 operator*(const Matrix3 &rhs) const noexcept;
    Vector operator*(const Vector &rhs) const noexcept;

    friend Vector;
    friend RotationMatrix;
    friend SkewMatrix;
    friend Matrix6;
    friend AdjungateMatrix;
    friend HomogenousMatrix;
};

Matrix3 operator*(double lhs, const Matrix3 &rhs) noexcept;

class RotationMatrix : public Matrix3 {
private:
    explicit RotationMatrix(Mat3  data) noexcept;
public:
    RotationMatrix(double z, double y, double x) noexcept;

    RotationMatrix operator*(const RotationMatrix &rhs) const noexcept;
    RotationMatrix inverse() const noexcept;

    friend AdjungateMatrix;
    friend HomogenousMatrix;
};

class SkewMatrix : public Matrix3 {
public:
    SkewMatrix(const Matrix3 & rhs);
    SkewMatrix(const Vector & vector) noexcept;
};

#endif //DAK_MATRIX3_H
