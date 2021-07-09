//
// Created by sba on 06.07.21.
//

#ifndef DAK_PLUECKER_H
#define DAK_PLUECKER_H

#include <eigen3/Eigen/Eigen>
#include "dual_number.h"

class Pluecker;
class SkewMatrix;
class Matrix3;
class Matrix6;
class AdjungateMatrix;

class PointVector;
class MomentVector;
class DirectionVector;
class HomogenousMatrix;

class Vector {
private:
    using Vec3 = Eigen::Matrix<double,3,1>;

    Vec3 data;

    explicit Vector(const Vec3 &data) noexcept;
public:
    Vector() = default;

    Vector(double a, double b, double c) noexcept;
    Vector(const SkewMatrix &skew) noexcept;

    Vector & operator+=(const Vector &rhs) noexcept;
    Vector operator+(const Vector &rhs) const noexcept;
    Vector operator-(const Vector &rhs) const noexcept;
    Vector operator*(double rhs) const noexcept;
    Vector operator/(double rhs) const;
    double operator*(const Vector &rhs) const noexcept; // dot
    Vector cross(const Vector &rhs) const noexcept;

    double norm() const noexcept;
    Vector normal() const;

    bool is_zero() const noexcept;

    friend Pluecker;
    friend PointVector;
    friend MomentVector;
    friend DirectionVector;
    friend Matrix3;
    friend SkewMatrix;
    friend HomogenousMatrix;
};

Vector operator*(double lhs, const Vector &rhs) noexcept;
Vector cross(const Vector &lhs, const Vector &rhs) noexcept;

struct ProjectionTrio;

class PointVector : public Vector {
public:
    explicit PointVector(const Vector &rhs) noexcept;
    PointVector(double a, double b, double c) noexcept;
};
class DirectionVector : public Vector {
public:
    explicit DirectionVector(const Vector &rhs);
    DirectionVector(double a, double b, double c);
};
class MomentVector : public Vector {
public:
    explicit MomentVector(const Vector &rhs) noexcept;
    MomentVector(double a, double b, double c) noexcept;
};

class Pluecker{
private:
    using Vec6 = Eigen::Matrix<double,6,1>;

    Vec6 data;
    bool is_valid;

    explicit Pluecker(Vec6 data) noexcept;
    Pluecker() = default;
public:
    explicit Pluecker(const PointVector &a, const PointVector &b) noexcept;
    explicit Pluecker(const DirectionVector &n, const PointVector &a) noexcept;
    explicit Pluecker(const DirectionVector &n, const MomentVector &m) noexcept;

    Vector n() const noexcept;
    Vector m() const noexcept;

    DualNumberAlgebra::DualNumber operator*(const Pluecker &rhs) const noexcept;

    Vector get_canonical_anchor() const noexcept;

    SkewMatrix get_direction_skew() const noexcept;

    AdjungateMatrix create_transform(DualNumberAlgebra::DualNumber angle) const noexcept;

    Pluecker operator-() const noexcept;
    Pluecker operator-(const Pluecker &rhs) const noexcept;

    Pluecker align() const;
    Pluecker normalize() const;

    ProjectionTrio project(const Pluecker &l) const noexcept;
    Vector point_project(const Pluecker &l) const noexcept;

    // not really line geometry
    Pluecker parallel_through_anchor(const Vector &new_anchor) const noexcept;

    friend Matrix6;
    friend AdjungateMatrix;
};

struct ProjectionTrio {
Pluecker p;
Pluecker r;
Pluecker o;
};

/*
class UnitPluecker {

};
 */

#endif //DAK_PLUECKER_H
