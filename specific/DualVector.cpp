//
// Created by sba on 29.06.20.
//

#include "DualVector.h"

DualNumberAlgebra::DualVector(const Vec3 &real) noexcept {
    this->_real = real;
    this->_dual << 0.0, 0.0, 0.0;
}

explicit DualNumberAlgebra::DualVector(const Vec3 &real, const Vec3 &dual) noexcept {
    this->_real = real;
    this->_dual = dual,
}

DualNumberAlgebra::DualVector &
DualNumberAlgebra::DualVector::operator=(const Vec3 &real) noexcept {
    this->_real = real;
    this->_dual << 0.0, 0.0, 0.0;
    return *this;
}

DualNumberAlgebra::DualVector
DualNumberAlgebra::DualVector::operator+(const DualVector &rhs) const noexcept {
    return DualNumberAlgebra::DualVector(this->_real + rhs._real, this->_dual + rhs._dual);
}

DualNumberAlgebra::DualVector
DualNumberAlgebra::DualVector::operator+(const Vec3 &real) const noexcept {
    return DualNumberAlgebra::DualVector(this->_real + real, this->_dual);
}

DualNumberAlgebra::DualVector
DualNumberAlgebra::DualVector::operator+() const noexcept {
    return DualNumberAlgebra::DualVector(this->_real, this->_dual);
}

DualNumberAlgebra::DualVector
DualNumberAlgebra::DualVector::operator-(const DualVector &rhs) const noexcept {
    return DualNumberAlgebra::DualVector(this->_real - rhs._real, this->_dual - rhs._dual);
}

DualNumberAlgebra::DualVector
DualNumberAlgebra::DualVector::operator-(const Vec3 &real) const noexcept {
    return DualNumberAlgebra::DualVector(this->_real - real, this->_dual);
}

DualNumberAlgebra::DualVector
DualNumberAlgebra::DualVector::operator-() const noexcept {
    return DualNumberAlgebra::DualVector(-this->_real, - this->_dual);
}

DualNumberAlgebra::DualVector
DualNumberAlgebra::DualVector::conjugate() const noexcept {
    return DualNumberAlgebra::DualVector(this->_real, - this->_dual);
}

double
DualNumberAlgebra::DualVector::norm() const noexcept {
    return this->_real.norm();
}

// projections = "getter"
Vec3
DualNumberAlgebra::DualVector::getReal() const noexcept {
    return this->_real;
}

Vec3
DualNumberAlgebra::DualVector::getDual() const noexcept {
    return this->_dual;
}

// mapping with derivatices

DualNumberAlgebra::DualNumber
DualNumberAlgebra::DualVector::operator*(const DualVector &rhs) const noexcept {
    return DualNumberAlgebra::DualNumber(this->_real.dot(rhs._real), this->_real.dot(rhs._dual) + rhs._real.dot(this->_dual));
}

DualNumberAlgebra::DualNumber
DualNumberAlgebra::DualVector::operator*(const Vec3 &real) const noexcept {
    return DualNumberAlgebra::DualNumber(this->_real.dot(real), this->_dual.dot(real));
}