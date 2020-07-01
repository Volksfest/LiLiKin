//
// Created by sba on 29.06.20.
//


#include <iomanip>

#include "RobotAlgebra.h"

RobotAlgebra::Matrix3::Matrix3() noexcept {
    data << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
}

RobotAlgebra::Matrix3::Matrix3(const RobotAlgebra::Vector3 &a, const RobotAlgebra::Vector3 &b,
                               const RobotAlgebra::Vector3 &c) noexcept : Matrix3(a.data, b.data, c.data)  {

}

RobotAlgebra::Matrix3::Matrix3(const Mat3 &mat) noexcept {
    data = mat;
}

RobotAlgebra::Matrix3::Matrix3(const Eigen::Matrix<double, 3, 1> &a, const Vec3 &b, const Vec3 &c) noexcept {
    data << a(0,0) , a(1,0) , a(2,0) ,
            b(0,0) , b(1,0) , b(2,0) ,
            c(0,0) , c(1,0) , c(2,0);
}

RobotAlgebra::Vector3::Vector3() noexcept {
    data << 0.0, 0.0, 0.0;
}

RobotAlgebra::Vector3::Vector3(const double &a, const double &b, const double &c) noexcept {
    data << a, b, c;
}

RobotAlgebra::Vector3::Vector3(const Vec3 &a) noexcept {
    data = a;
}

RobotAlgebra::Matrix3
RobotAlgebra::Vector3::cross_mat() const noexcept {
    Mat3 m;
    m << 0.0, -data(2,0), data(1,0),
            data(2,0), 0.0, -data(0,0),
            -data(1,0), data(0,0), 0.0;
    return Matrix3(m);
}

RobotAlgebra::Adjungate::Adjungate(const RobotAlgebra::Vector3 &n, const RobotAlgebra::Vector3 &m) noexcept {
    data.topLeftCorner(3, 3)     = n.cross_mat().data;
    data.topRightCorner(3, 3)     = Mat3::Zero(3, 3);
    data.bottomLeftCorner(3, 3)     = m.cross_mat().data;
    data.bottomRightCorner(3, 3)     = n.cross_mat().data;
}

RobotAlgebra::Pluecker::Pluecker(const RobotAlgebra::Vector3 &n, const RobotAlgebra::Vector3 &m) noexcept {
    this->n = n;
    this->m = m;
}

RobotAlgebra::Pluecker
RobotAlgebra::Pluecker::fromPoints(const RobotAlgebra::Vector3 &a, const RobotAlgebra::Vector3 &b) {
    Vec3 n = b.data - a.data;
    n = n / n.norm();
    Vec3 m = a.data.cross(n);
    return Pluecker(Vector3(n),Vector3(m));
}

std::ostream &std::operator<<(std::ostream &stream, RobotAlgebra::Matrix3 const &d) {
    return stream << d.data << std::endl;
}

std::ostream &std::operator<<(std::ostream &stream, RobotAlgebra::Vector3 const &d) {
    return stream << d.data << std::endl;
}

std::ostream &std::operator<<(std::ostream &stream, RobotAlgebra::Pluecker const &d) {
    return stream << d.n << std::endl << d.m << std::endl;
}

