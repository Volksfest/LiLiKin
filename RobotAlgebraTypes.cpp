#include <utility>

//
// Created by sba on 17.09.20.
//

#include "RobotAlgebraTypes.h"

RobotAlgebra::SkewMatrix RobotAlgebra::SkewMatrix::from(const RobotAlgebra::Vector &vector) noexcept {
    SkewMatrix m;
    m.data << 0, -vector.data(2, 0), vector.data(1, 0),
            vector.data(2, 0), 0, -vector.data(0, 0),
            -vector.data(1, 0), vector.data(0, 0), 0;
    return m;
}

RobotAlgebra::RotationMatrix::RotationMatrix(double z, double y, double x) noexcept {
    Eigen::AngleAxisd rollAngle(z, Eigen::Vector3d::UnitZ());
    Eigen::AngleAxisd yawAngle(y, Eigen::Vector3d::UnitY());
    Eigen::AngleAxisd pitchAngle(x, Eigen::Vector3d::UnitX());

    Eigen::Quaternion<double> q = rollAngle * yawAngle * pitchAngle;

    this->data = q.matrix();
}

RobotAlgebra::RotationMatrix
RobotAlgebra::RotationMatrix::inverse() const noexcept {

    return RotationMatrix(this->data.transpose());
}

RobotAlgebra::GenericMatrix3::GenericMatrix3(Eigen::Matrix<double,3,3>  data) noexcept : data(std::move(data)) { }

RobotAlgebra::GenericMatrix3
RobotAlgebra::GenericMatrix3::operator*(
        double rhs) const noexcept {

    return GenericMatrix3(this->data * rhs);
}

RobotAlgebra::GenericMatrix3
RobotAlgebra::GenericMatrix3::operator*(
        const RobotAlgebra::GenericMatrix3 & rhs) const noexcept {

    return GenericMatrix3(this->data * rhs.data);
}

RobotAlgebra::Vector
RobotAlgebra::GenericMatrix3::operator*(const Vector &rhs) const noexcept {
    return Vector(this->data * rhs.data);
}

RobotAlgebra::GenericMatrix3 RobotAlgebra::operator*(double lhs, const RobotAlgebra::GenericMatrix3 &rhs) noexcept {
    return rhs * lhs;
}

RobotAlgebra::Vector::Vector(Eigen::Matrix<double, 3, 1> data) noexcept : data(std::move(data)) { }

RobotAlgebra::Vector::Vector(double a, double b, double c) noexcept {
    this->data << a, b, c;
}

RobotAlgebra::Vector
RobotAlgebra::Vector::operator+(
        const RobotAlgebra::Vector &rhs) const noexcept {

    return Vector(this->data + rhs.data);
}

RobotAlgebra::Vector
RobotAlgebra::Vector::operator-(
        const RobotAlgebra::Vector &rhs) const noexcept {

    return Vector(this->data - rhs.data);
}

double
RobotAlgebra::Vector::operator*(
        const RobotAlgebra::Vector &rhs) const noexcept {

    return this->data.dot(rhs.data);
}

RobotAlgebra::Vector operator*(double lhs, const RobotAlgebra::Vector &rhs) noexcept {
    return rhs * lhs;
}

RobotAlgebra::Vector
RobotAlgebra::Vector::cross(
        const RobotAlgebra::Vector &rhs) const noexcept {

    return Vector(this->data.cross(rhs.data));
}

double
RobotAlgebra::Vector::norm() const noexcept {

    return this->data.norm();
}

RobotAlgebra::Vector
RobotAlgebra::Vector::operator*(
        double rhs) const noexcept {

    return Vector(this->data * rhs);
}

RobotAlgebra::Vector
RobotAlgebra::Vector::operator/(
        double rhs) const noexcept {

    return Vector(this->data / rhs);
}

RobotAlgebra::Vector
RobotAlgebra::Vector::from(
        const SkewMatrix &skew) noexcept {

    Vector n;
    n.data << skew.data(2,1), skew.data(0,2), skew.data(1,0);
    return n;
}

RobotAlgebra::Pluecker::Pluecker(Eigen::Matrix<double, 6, 1> data) noexcept : data(std::move(data)) { }

RobotAlgebra::Pluecker
RobotAlgebra::Pluecker::fromPoints(
        const Vector &a,
        const Vector &b) noexcept {

    Vector n = b - a;
    n = n / n.norm();
    Vector m = a.cross(n);
    Pluecker p;
    p.data << n.data, m.data;
    return p;
}

RobotAlgebra::Pluecker
RobotAlgebra::Pluecker::fromDirection(
        const Vector &n,
        const Vector &a) noexcept {

    Vector m = a.cross(n);
    Pluecker p;
    p.data << n.data, m.data;
    return p;
}

RobotAlgebra::Vector
RobotAlgebra::Pluecker::v() const noexcept {
    return Vector(this->data.head(3));
}

RobotAlgebra::Vector
RobotAlgebra::Pluecker::w() const noexcept {
    return Vector(this->data.tail(3));
}

RobotAlgebra::AdjungateMatrix
RobotAlgebra::Pluecker::transform(
        DualNumberAlgebra::DualNumber<double> angle) const noexcept {

    GenericMatrix6 orthoterm = GenericMatrix6::from(*this);
    GenericMatrix6 uniterm = - orthoterm * orthoterm; //eq -i^2
    GenericMatrix6 nullterm(Eigen::Matrix<double,6,6>::Identity(6,6) - uniterm.data);

    return AdjungateMatrix((GenericMatrix6::from(cos(angle)) * uniterm +
           GenericMatrix6::from(sin(angle)) * orthoterm +
           nullterm).data);
}

RobotAlgebra::Pluecker
RobotAlgebra::Pluecker::normalize() const noexcept {

    double norm = this->data.head(3).norm();
    // TODO epsilon
    if(norm == 0) {
        // TODO ? Throw exception, thats not a line...
        return Pluecker();
    }

    return Pluecker(this->data / norm); // TODO all?
}

RobotAlgebra::Pluecker
RobotAlgebra::Pluecker::align() const noexcept {
    Pluecker a;

    Eigen::Matrix<double,3,1> n,m;
    n = this->data.head(3);
    m = this->data.tail(3);

    a.data.head(3) = n;
    a.data.tail(3) = m - (n.dot(m) / n.dot(n)) * n;
    return a;
}

void
RobotAlgebra::Pluecker::project(
        const RobotAlgebra::Pluecker &n,
        RobotAlgebra::Pluecker &p,
        RobotAlgebra::Pluecker &r,
        RobotAlgebra::Pluecker &o) const noexcept {

    GenericMatrix6 adj = GenericMatrix6::from(n);

    o = Pluecker(adj.data * this->data);
    r = Pluecker(-adj.data * o.data);
    p = Pluecker(this->data - r.data); // TODO either check unity or multiply with norm
    if (o.v().norm() < 1) {
        auto anchor_n = n.v().cross(n.w());
        auto c = anchor_n - this->v().cross(this->w());
        o = Pluecker::fromDirection(n.v().cross(c).cross(n.v()), anchor_n).normalize().align();
    }
}

RobotAlgebra::RotationMatrix
RobotAlgebra::HomogenousMatrix::R() const noexcept {

    return RotationMatrix(this->data.topLeftCorner(3,3));
}

RobotAlgebra::Vector
RobotAlgebra::HomogenousMatrix::p() const noexcept {

    return Vector(this->data.topRightCorner(3,1));
}

RobotAlgebra::HomogenousMatrix
RobotAlgebra::HomogenousMatrix::from(
        const RobotAlgebra::RotationMatrix &rot,
        const RobotAlgebra::Vector &vector) noexcept {

    HomogenousMatrix hom;
    hom.data = Eigen::Matrix<double,4,4>::Zero(4,4);

    hom.data.topLeftCorner(3,3) = rot.data;
    hom.data.topRightCorner(3,1) = vector.data;
    hom.data(3,3) = 1.0;

    return hom;
}

RobotAlgebra::HomogenousMatrix
RobotAlgebra::HomogenousMatrix::from(
        const RobotAlgebra::AdjungateMatrix & adj) noexcept {


    RotationMatrix rot = adj.R();
    GenericMatrix3 trans_rot = adj.pxR();
    //trans_rot.data = adj.data.bottomLeftCorner(3,3);
    SkewMatrix sk( (trans_rot * rot.inverse()).data);

    return HomogenousMatrix::from(
            rot,
            Vector::from(sk));
}

RobotAlgebra::GenericMatrix6::GenericMatrix6(Eigen::Matrix<double, 6, 6> data) noexcept :data(std::move(data)){ }

RobotAlgebra::GenericMatrix6
RobotAlgebra::GenericMatrix6::from(
        const RobotAlgebra::GenericMatrix3 &real,
        const RobotAlgebra::GenericMatrix3 &dual) noexcept {

    GenericMatrix6 mat;
    mat.data.topLeftCorner(3, 3)     = real.data;
    mat.data.topRightCorner(3, 3)     = Eigen::Matrix<double,3,3>::Zero(3, 3);
    mat.data.bottomLeftCorner(3, 3)     = dual.data;
    mat.data.bottomRightCorner(3, 3)     = real.data;
    return mat;
}

RobotAlgebra::GenericMatrix6
RobotAlgebra::GenericMatrix6::from(
        const DualNumberAlgebra::DualNumber<double> &dn) noexcept {

    GenericMatrix3 i3(Eigen::Matrix<double,3,3>::Identity(3,3));
    return RobotAlgebra::GenericMatrix6::from(
            dn.getReal() * i3,
            dn.getDual() * i3);
}

RobotAlgebra::GenericMatrix6
RobotAlgebra::GenericMatrix6::from(
        const Pluecker &pluecker) noexcept {

    auto v = pluecker.v();
    auto w = pluecker.w();
    return GenericMatrix6::from(
            SkewMatrix::from(v),
            SkewMatrix::from(w));
}

RobotAlgebra::GenericMatrix6
RobotAlgebra::GenericMatrix6::operator+(
        const RobotAlgebra::GenericMatrix6 &rhs) const noexcept {

    return RobotAlgebra::GenericMatrix6(this->data + rhs.data);
}

RobotAlgebra::GenericMatrix6
RobotAlgebra::GenericMatrix6::operator-(
        const RobotAlgebra::GenericMatrix6 &rhs) const noexcept {

    return RobotAlgebra::GenericMatrix6(this->data - rhs.data);
}

RobotAlgebra::GenericMatrix6
RobotAlgebra::GenericMatrix6::operator*(
        const RobotAlgebra::GenericMatrix6 &rhs) const noexcept {

    return RobotAlgebra::GenericMatrix6(this->data * rhs.data);
}

RobotAlgebra::GenericMatrix6
RobotAlgebra::GenericMatrix6::operator-() const noexcept {

    return RobotAlgebra::GenericMatrix6(-this->data);
}

RobotAlgebra::Pluecker
RobotAlgebra::GenericMatrix6::operator*(
        const RobotAlgebra::Pluecker &rhs) const noexcept {

    return RobotAlgebra::Pluecker(this->data * rhs.data);
}

RobotAlgebra::AdjungateMatrix
RobotAlgebra::AdjungateMatrix::from(
        const RobotAlgebra::RotationMatrix &real,
        const RobotAlgebra::GenericMatrix3 &dual) noexcept {

    AdjungateMatrix mat;
    mat.data.topLeftCorner(3, 3)     = real.data;
    mat.data.topRightCorner(3, 3)     = Eigen::Matrix<double,3,3>::Zero(3, 3);
    mat.data.bottomLeftCorner(3, 3)     = dual.data;
    mat.data.bottomRightCorner(3, 3)     = real.data;
    return mat;
}

RobotAlgebra::RotationMatrix
RobotAlgebra::AdjungateMatrix::R() const noexcept {

    return RotationMatrix(this->data.topLeftCorner(3,3));
}

RobotAlgebra::GenericMatrix3
RobotAlgebra::AdjungateMatrix::pxR() const noexcept {

    return GenericMatrix3(this->data.bottomLeftCorner(3,3));
}

RobotAlgebra::AdjungateMatrix
RobotAlgebra::AdjungateMatrix::from(
        const HomogenousMatrix &hom) noexcept {

    return from(
            hom.R(),
            SkewMatrix::from(hom.p()) * hom.R());
}
