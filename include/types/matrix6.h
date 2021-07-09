//
// Created by sba on 06.07.21.
//

#ifndef DAK_MATRIX6_H
#define DAK_MATRIX6_H

#include <eigen3/Eigen/Eigen>

#include "dual_number.h"

class Matrix3;
class Pluecker;
class RotationMatrix;
class HomogenousMatrix;

class Vector;
class AdjungateMatrix;

class Matrix6 {
private:
    using Mat6 = Eigen::Matrix<double, 6, 6>;

    Eigen::Matrix<double,6,6> data;

    explicit Matrix6(Mat6 data) noexcept;
public:
    Matrix6(const Matrix3 &real, const Matrix3 &dual) noexcept;
    explicit Matrix6(const DualNumberAlgebra::DualNumber &dn) noexcept;
    explicit Matrix6(const Pluecker &pl) noexcept;

    Matrix6 operator-() const noexcept;
    Matrix6 operator+(const Matrix6 &rhs) const noexcept;
    Matrix6 operator-(const Matrix6 &rhs) const noexcept;
    Matrix6 operator*(const Matrix6 &rhs) const noexcept;
    Pluecker operator*(const Pluecker &rhs) const noexcept;

    friend AdjungateMatrix;
};

class AdjungateMatrix : public Matrix6 {
private:
    explicit AdjungateMatrix(const Mat6 &mat) noexcept;
    explicit AdjungateMatrix(const Matrix6 &mat) noexcept;
public:
    // To constructor
    AdjungateMatrix(const RotationMatrix &real, const Matrix3 &dual) noexcept;
    explicit AdjungateMatrix(const HomogenousMatrix &pl) noexcept;


    AdjungateMatrix operator*(const AdjungateMatrix &rhs) const noexcept;
    Pluecker operator*(const Pluecker &rhs) const noexcept;

    AdjungateMatrix inverse() const noexcept;

    RotationMatrix R() const noexcept;
    Matrix3 pxR() const noexcept;
    Vector p() const noexcept;

    friend Pluecker;
};

#endif //DAK_MATRIX6_H
