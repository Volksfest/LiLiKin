//
// Created by sba on 29.06.20.
//

#include "DualNumber.h"

DualNumberAlgebra::DualNumber::DualNumber(const double &real, const double &dual) noexcept {
    this->_real = real;
    this->_dual = dual;
}
/*
DualNumberAlgebra::DualNumber &
DualNumberAlgebra::DualNumber::operator=(const DualNumberAlgebra::DualNumber &rhs) noexcept {
    this->_real = rhs._real;
    this->_dual = rhs._dual;
    return *this;
}*/

DualNumberAlgebra::DualNumber &
DualNumberAlgebra::DualNumber::operator=(const double &real) noexcept {
    this->_real = real;
    this->_dual = 0;
    return *this;
}

DualNumberAlgebra::DualNumber
DualNumberAlgebra::DualNumber::operator+(const DualNumberAlgebra::DualNumber &rhs) const noexcept {
    return DualNumber(this->_real + rhs._real, this->_dual + rhs._dual);
}

DualNumberAlgebra::DualNumber
DualNumberAlgebra::DualNumber::operator+(const double &real) const noexcept {
    return DualNumber(this->_real + real, this->_dual);
}

DualNumberAlgebra::DualNumber
DualNumberAlgebra::DualNumber::operator+() const noexcept {
    return DualNumber(this->_real, this->_dual);
}

DualNumberAlgebra::DualNumber
DualNumberAlgebra::DualNumber::operator*(const DualNumberAlgebra::DualNumber &rhs) const noexcept {
    return DualNumber(this->_real * rhs._real, this->_dual * rhs._real + this->_real * rhs._dual);
}

DualNumberAlgebra::DualNumber
DualNumberAlgebra::DualNumber::operator*(const double &real) const noexcept {
    return DualNumber(this->_real * real, this->_dual);
}

DualNumberAlgebra::DualNumber
DualNumberAlgebra::DualNumber::operator-(const DualNumberAlgebra::DualNumber &rhs) const noexcept {
    return DualNumber(this->_real - rhs._real, this->_dual - rhs._dual);
}

DualNumberAlgebra::DualNumber
DualNumberAlgebra::DualNumber::operator-(const double &real) const noexcept {
    return DualNumber(this->_real - real, this->_dual);
}

DualNumberAlgebra::DualNumber
DualNumberAlgebra::DualNumber::operator-() const noexcept {
    return DualNumber(-this->_real, -this->_dual);
}

DualNumberAlgebra::DualNumber
DualNumberAlgebra::DualNumber::conjugate() const noexcept {
    return DualNumber(this->_real, -this->_dual);
}

DualNumberAlgebra::DualNumber
DualNumberAlgebra::DualNumber::inverse() const {
    if (this->_real == 0) {
        throw std::logic_error("Cannot invert a dual number with real part equals to zero");
    }
    return DualNumber(1 / this->_real, -this->_dual / this->_real / this->_real);
}

double
DualNumberAlgebra::DualNumber::norm() const noexcept {
    return this->_real >= 0 ? this->_real : -this->_real;
}

DualNumberAlgebra::DualNumber
DualNumberAlgebra::DualNumber::operator/(const DualNumberAlgebra::DualNumber &rhs) const {
    if (this->_real == 0) {
        throw std::logic_error("Cannot divide by a dual number with real part equals to zero");
    }
    return DualNumber(this->_real / rhs._real,
                      this->_dual / rhs._real - this->_real * rhs._dual / rhs._real / rhs._real);
}

double
DualNumberAlgebra::DualNumber::getReal() const noexcept {
    return this->_real;
}

double
DualNumberAlgebra::DualNumber::getDual() const noexcept {
    return this->_dual;
}

DualNumberAlgebra::DualNumber
DualNumberAlgebra::DualNumber::sin(const DualNumberAlgebra::DualNumber &phi) const noexcept {
    const double &r = phi.getReal();
    const double &d = phi.getDual();
    return DualNumber(std::sin(r), d * std::cos(r));
}

DualNumberAlgebra::DualNumber
DualNumberAlgebra::DualNumber::asin(const DualNumberAlgebra::DualNumber &w) const noexcept {
    const double &r = w.getReal();
    const double nr = std::asin(r);
    const double &d = w.getDual();
    return DualNumber(nr, d / std::cos(nr));
}

DualNumberAlgebra::DualNumber
DualNumberAlgebra::DualNumber::cos(const DualNumberAlgebra::DualNumber &phi) const noexcept {
    const double &r = phi.getReal();
    const double &d = phi.getDual();
    return DualNumber(std::cos(r), -d * std::sin(r));
}

DualNumberAlgebra::DualNumber
DualNumberAlgebra::DualNumber::acos(const DualNumberAlgebra::DualNumber &w) const noexcept {
    const double &r = w.getReal();
    const double nr = std::acos(r);
    const double &d = w.getDual();
    return DualNumber(nr, -d / std::sin(nr));
}

DualNumberAlgebra::DualNumber
DualNumberAlgebra::DualNumber::tan(const DualNumberAlgebra::DualNumber &phi) const noexcept {
    const double &r = phi.getReal();
    const double &d = phi.getDual();
    const double cr = std::cos(r);
    return DualNumber(std::tan(r), d / cr / cr);
}

DualNumberAlgebra::DualNumber
DualNumberAlgebra::DualNumber::atan(const DualNumberAlgebra::DualNumber &w) const noexcept {
    const double &r = w.getReal();
    const double nr = std::atan(r);
    const double &d = w.getDual();
    const double cnr = std::cos(nr);
    return DualNumber(nr, -d * cnr * cnr);
}

DualNumberAlgebra::DualNumber
DualNumberAlgebra::DualNumber::atan2(const DualNumberAlgebra::DualNumber &y,
                                                                   const DualNumberAlgebra::DualNumber &x) const noexcept {
    const double &xr = x.getReal();
    const double &yr = y.getReal();
    const double &xd = x.getDual();
    const double &yd = y.getDual();
    return DualNumber(std::atan2(yr, xr), (xr*yd - xd * yr) / ( xr * xr + yr * yr));
}

DualNumberAlgebra::DualNumber
DualNumberAlgebra::operator+(const double &lhs, const DualNumberAlgebra::DualNumber &rhs) noexcept {
    return rhs + lhs;
}

DualNumberAlgebra::DualNumber
DualNumberAlgebra::operator ""_s(long double dual) noexcept {
    return DualNumber(0.0, dual);
}

std::ostream &std::operator<<(std::ostream &stream, DualNumberAlgebra::DualNumber const &d) {
    double dual = d.getDual();
    return (stream << "(" << d.getReal() << (dual < 0 ? "-" : "+") << (dual < 0 ? -dual : dual) << "s)");
}
