//
// Created by sba on 07.07.21.
//

#include "base/vector.h"
#include "base/matrix3.h"

Vector::Vector(const Vec3 &data) noexcept {
    this->data = data;
}

Vector::Vector(double a, double b, double c) noexcept {
    this->data << a,b,c;
}

Vector::Vector(const SkewMatrix &skew) noexcept {
    this->data << skew.get()(2,1), skew.get()(0,2), skew.get()(1,0);
}

Vector &
Vector::operator+=(const Vector &rhs) noexcept {
    this->data += rhs.data;
    return *this;
}

Vector
Vector::operator+(const Vector &rhs) const noexcept {
    return Vector(this->data + rhs.data);
}

Vector
Vector::operator-(const Vector &rhs) const noexcept {
    return Vector(this->data - rhs.data);
}

Vector
Vector::operator*(double rhs) const noexcept {
    return Vector(this->data *rhs);
}

Vector
operator*(double lhs, const Vector &rhs) noexcept {
    return rhs * lhs;
}

Vector
Vector::operator/(double rhs) const {
    if (rhs == 0.0) {
        throw std::domain_error("Division by zero");
    }
    return Vector(this->data / rhs);
}

double
Vector::operator*(const Vector &rhs) const noexcept {
    return this->data.dot(rhs.data);
}

Vector
Vector::cross(const Vector &rhs) const noexcept {
    return Vector(this->data.cross(rhs.data));
}

double
Vector::norm() const noexcept {
    return this->data.norm();
}

Vector
Vector::normal() const {
    return *this / this->norm();
}

bool
Vector::is_zero() const noexcept {
    return this->norm() < 0.00001;
}

Vector cross(const Vector &lhs, const Vector &rhs) noexcept {
    return lhs.cross(rhs);
}

const Vector::Vec3 &
Vector::get() const noexcept {
    return data;
}

bool
Vector::operator==(const Vector &rhs) const noexcept {
    return this->data.isApprox(rhs.data);
}

PointVector::PointVector(const Vector &rhs) noexcept : Vector(rhs) {}

PointVector::PointVector(double a, double b, double c) noexcept : Vector(a,b,c) {}

DirectionVector::DirectionVector(const Vector &rhs) : Vector(rhs) {
    if (this->is_zero()) {
        throw std::invalid_argument("Given vector is the null vector");
    }
}

DirectionVector::DirectionVector(double a, double b, double c) : Vector(a,b,c) {
    if (this->is_zero()) {
        throw std::invalid_argument("Given vector is the null vector");
    }
}

UnitDirectionVector
DirectionVector::normal() const {
    return UnitDirectionVector(*this / this->norm());
}

UnitDirectionVector::UnitDirectionVector(const Vector &rhs) : DirectionVector(rhs) {
    if (abs(1.0 - this->norm()) > 0.0001) { // TODO epsilon
        throw std::invalid_argument("Given vector is not an unit vector");
    }
}

UnitDirectionVector::UnitDirectionVector(double a, double b, double c) : DirectionVector(a,b,c) {
    if (abs(1.0 - this->norm()) > 0.0001) { // TODO epsilon
        throw std::invalid_argument("Given vector is not an unit vector");
    }
}

UnitDirectionVector
UnitDirectionVector::normal() const {
    return *this;
}


MomentVector::MomentVector(const Vector &rhs) noexcept : Vector(rhs) {}

MomentVector::MomentVector(double a, double b, double c) noexcept : Vector(a,b,c) {}