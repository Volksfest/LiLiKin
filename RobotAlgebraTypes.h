//
// Created by sba on 17.09.20.
//

#ifndef DAK_ROBOTALGEBRATYPES_H
#define DAK_ROBOTALGEBRATYPES_H

#include <eigen3/Eigen/Eigen>
#include "DualNumber.h"

namespace RobotAlgebra {

    class SkewMatrix;
    class Vector {
    public:
        Vector() = default;
        explicit Vector(Eigen::Matrix<double,3,1>  data) noexcept;
        Vector(double a, double b, double c) noexcept;
        Eigen::Matrix<double,3,1> data;
        Vector operator+(const Vector &rhs) const noexcept;
        Vector operator-(const Vector &rhs) const noexcept;
        double operator*(const Vector &rhs) const noexcept;
        Vector cross(const Vector &rhs) const noexcept;
        double norm() const noexcept;
        Vector operator*(double rhs) const noexcept;
        Vector operator/(double rhs) const noexcept;

        static Vector from(const SkewMatrix &skew) noexcept;
    };

    Vector operator*(double lhs, const Vector &rhs) noexcept;

    class AdjungateMatrix;
    class Pluecker{
    public:
        Pluecker() = default;
        explicit Pluecker(Eigen::Matrix<double,6,1>  data) noexcept;
        static Pluecker fromPoints(const Vector &a, const Vector &b) noexcept;
        static Pluecker fromDirection(const Vector &n, const Vector &a) noexcept;
        Eigen::Matrix<double,6,1> data;

        Vector v() const noexcept;
        Vector w() const noexcept;

        AdjungateMatrix transform(DualNumberAlgebra::DualNumber<double> angle) const noexcept;

        Pluecker align() const noexcept;
        Pluecker normalize() const noexcept;
        void project(const Pluecker &n, Pluecker &p, Pluecker &r, Pluecker &o) const noexcept;
    };

    class GenericMatrix3 {
    public:
        GenericMatrix3() = default;
        explicit GenericMatrix3(Eigen::Matrix<double,3,3>  data) noexcept;

        GenericMatrix3 operator*(double rhs) const noexcept;
        GenericMatrix3 operator*(const GenericMatrix3 &rhs) const noexcept;
        Vector operator*(const Vector &rhs) const noexcept;

        Eigen::Matrix<double,3,3> data;
    };

    GenericMatrix3 operator*(double lhs, const GenericMatrix3 &rhs) noexcept;


    class RotationMatrix : public GenericMatrix3 {
        // TODO change visibility
        using GenericMatrix3::GenericMatrix3;
    public:
        RotationMatrix(double z, double y, double z2) noexcept;
        RotationMatrix inverse() const noexcept;
    };

    class SkewMatrix : public GenericMatrix3 {
        using GenericMatrix3::GenericMatrix3;
    public:
        static SkewMatrix from(const Vector & vector) noexcept;
    };


    class AdjungateMatrix;

    class HomogenousMatrix {
    public:
        Eigen::Matrix<double,4,4> data;

        static HomogenousMatrix from(const RotationMatrix & rot, const Vector & vector) noexcept;
        static HomogenousMatrix from(const AdjungateMatrix &adj) noexcept;

        RotationMatrix R() const noexcept;
        Vector p() const noexcept;
    };

    class GenericMatrix6 {
    public:
        GenericMatrix6() = default;
        explicit GenericMatrix6(Eigen::Matrix<double, 6, 6> data) noexcept;
        Eigen::Matrix<double,6,6> data;
        static GenericMatrix6 from(const GenericMatrix3 &real, const GenericMatrix3 &dual) noexcept;
        static GenericMatrix6 from(const DualNumberAlgebra::DualNumber<double> &dn) noexcept;

        static GenericMatrix6 from(const Pluecker &pl) noexcept;

        GenericMatrix6 operator-() const noexcept;
        GenericMatrix6 operator+(const GenericMatrix6 &rhs) const noexcept;
        GenericMatrix6 operator-(const GenericMatrix6 &rhs) const noexcept;
        GenericMatrix6 operator*(const GenericMatrix6 &rhs) const noexcept;
        Pluecker operator*(const Pluecker &rhs) const noexcept;
    };

    class AdjungateMatrix : public GenericMatrix6 {
        using GenericMatrix6::GenericMatrix6;
    public:
        static AdjungateMatrix from(const RotationMatrix &real, const GenericMatrix3 &dual) noexcept;
        static AdjungateMatrix from(const HomogenousMatrix &pl) noexcept;

        RotationMatrix R() const noexcept;
        GenericMatrix3 pxR() const noexcept;
    };
}


#endif //DAK_ROBOTALGEBRATYPES_H
