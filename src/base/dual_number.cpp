//
// Created by sba on 06.07.21.
//

#include "base/dual_number.h"
#include <iostream>

namespace DualNumberAlgebra {

// function mapping so other functions can be mapped
    double (*real_sin)(double) = std::sin;

    double (*real_asin)(double) = std::asin;

    double (*real_cos)(double) = std::cos;

    double (*real_acos)(double) = std::acos;

    double (*real_tan)(double) = std::tan;

    double (*real_atan)(double) = std::atan;

    double (*real_atan2)(double, double) = std::atan2;

    double (*real_sqrt)(double) = std::sqrt;

    DualNumber &
    DualNumber::operator=(const double &real) noexcept {
        this->_real = real;
        this->_dual = 0;
        return *this;
    }

    DualNumber
    DualNumber::operator+(const DualNumber &rhs) const noexcept {
        return DualNumber(this->_real + rhs._real, this->_dual + rhs._dual);
    }

    DualNumber
    DualNumber::operator+() const noexcept {
        return DualNumber(this->_real, this->_dual);
    }

    DualNumber &
    DualNumber::operator+=(const DualNumber &rhs) noexcept {
        this->_real += rhs._real;
        this->_dual += rhs._dual;
        return *this;
    }

    DualNumber
    DualNumber::operator*(const DualNumber &rhs) const noexcept {
        return DualNumber(this->_real * rhs._real, this->_dual * rhs._real + this->_real * rhs._dual);
    }

    DualNumber &
    DualNumber::operator*=(const DualNumber &rhs) noexcept {
        this->_real *= rhs._real;
        this->_dual *= rhs._dual;
        return *this;
    }


    DualNumber
    DualNumber::operator-(const DualNumber &rhs) const noexcept {
        return DualNumber(this->_real - rhs._real, this->_dual - rhs._dual);
    }

    DualNumber
    DualNumber::operator-() const noexcept {
        return DualNumber(-this->_real, -this->_dual);
    }

    DualNumber &
    DualNumber::operator-=(const DualNumber &rhs) noexcept {
        this->_real -= rhs._real;
        this->_dual -= rhs._dual;
        return *this;
    }

    DualNumber
    DualNumber::operator/(const DualNumber &rhs) const {
        if (rhs._real == 0) {
            throw std::logic_error("Cannot divide by (a real part equal to) zero");
        }
        return DualNumber(this->_real / rhs._real,
                          this->_dual / rhs._real - this->_real * rhs._dual / rhs._real / rhs._real);
    }

    DualNumber &
    DualNumber::operator/=(const DualNumber &rhs) {
        if (rhs._real == 0) {
            throw std::logic_error("Cannot divide by (a real part equal to) zero");
        }
        this->_real /= rhs._real;
        this->_dual /= rhs._real;
        this->_dual -= this->_real * rhs._dual / rhs._real / rhs._real;
        return *this;
    }

    DualNumber
    DualNumber::conjugate() const noexcept {
        return DualNumber(this->_real, -this->_dual);
    }

    DualNumber
    DualNumber::inverse() const {
        if (this->_real == 0) {
            throw std::logic_error("Cannot invert (a dual number with real part equals to) zero");
        }
        return DualNumber(1 / this->_real, -this->_dual / this->_real / this->_real);
    }

    double
    DualNumber::norm() const noexcept {
        return this->_real >= 0 ? this->_real : -this->_real;
    }

    double
    DualNumber::norm_square() const noexcept {
        return this->_real * this->_real;
    }

    double
    DualNumber::real() const noexcept {
        return this->_real;
    }

    double
    DualNumber::dual() const noexcept {
        return this->_dual;
    }

    bool
    DualNumber::is_zero() const noexcept {
        return this->_real == 0 && this->_dual == 0; //TODO epsilon
    }

    DualNumber
    sqrt(const DualNumber &x) noexcept {
        auto root = real_sqrt(x.real());
        return DualNumber(
                root,
                0.5 * x.dual() / root);
    }

    DualNumber
    sin(const DualNumber &phi) noexcept {
        return DualNumber(
                real_sin(phi.real()),
                phi.dual() * real_cos(phi.real()));
    }

    DualNumber
    asin(const DualNumber &w) noexcept {
        return DualNumber(
                real_asin(w.real()),
                w.dual() / real_sqrt(1 - w.real() * w.real()));
    }

    DualNumber
    cos(const DualNumber &phi) noexcept {
        return DualNumber(
                real_cos(phi.real()),
                -phi.dual() * real_sin(phi.real()));
    }

    DualNumber
    acos(const DualNumber &w) noexcept {
        return DualNumber(
                real_acos(w.real()),
                -w.dual() / real_sqrt(1 - w.real() * w.real()));
    }

    DualNumber
    tan(const DualNumber &phi) noexcept {
        const double cr = real_cos(phi.real());

        return DualNumber(
                real_tan(phi.real()),
                phi.dual() / cr / cr);
    }

    DualNumber
    atan(const DualNumber &w) noexcept {
        return DualNumber(
                real_atan(w.real()),
                w.dual() / (1 + w.real() * w.real()));
    }

    DualNumber
    atan2(const DualNumber &y, const DualNumber &x) noexcept {
        const double &xr = x.real();
        const double &yr = y.real();
        const double &xd = x.dual();
        const double &yd = y.dual();
        return DualNumber(
                real_atan2(yr, xr),
                (xr * yd - xd * yr) /
                (xr * xr + yr * yr));
    }

    std::ostream &operator<<(std::ostream &stream, DualNumber const &d) {
        double dual = d.dual();
        return (stream << d.real() << (dual < 0 ? "-" : "+") << std::abs(dual) << "ϵ");
    }

    std::vector<DualNumber>
    solve_trigonometric_equation(const DualNumber &cos_factor,const DualNumber &sin_factor, const DualNumber &offset) {
        DualNumber dd = cos_factor * cos_factor +
                        sin_factor * sin_factor -
                        offset * offset;

        if(abs(dd.real()) > 0.000000001) {
            if (dd.real() < 0) {
                throw std::domain_error("No solution possible");
            }
        }

        DualNumber pre = atan2(sin_factor, cos_factor);

        // TODO epsilon
        if(dd.real() < 0.000000001) {
            // there is some problems with calculation d = sqrt(dd) if dd has a zero real part
            // but luckily d is not necessary if the real part is zero
            return {pre};
        } else {
            DualNumber d = DualNumberAlgebra::sqrt(dd);
            DualNumber rad = atan2(d, offset);
            return {pre + rad, pre - rad};
        }
    }

    DualNumber
    operator+(double lhs, DualNumber rhs) noexcept {
        return rhs + lhs;
    }

    DualNumber
    operator*(double lhs, DualNumber rhs) noexcept {
        return rhs * lhs;
    }

    DualNumber
    operator-(double lhs, DualNumber rhs) noexcept {
        return DualNumber(lhs) - rhs;
    }

    DualNumber
    operator/(double lhs, DualNumber rhs) {
        return DualNumber(lhs) / rhs;
    }

    bool
    operator==(const DualNumber &lhs, const DualNumber &rhs) noexcept {
        return (lhs.real() == rhs.real()) && (lhs.dual() == rhs.dual());
    }

    namespace literals{
        DualNumber
        operator ""_s(long double dual) noexcept {
            return DualNumber(0.0, dual);
        }

        DualNumber
        operator ""_s(unsigned long long dual) noexcept {
            return DualNumber(0.0, static_cast<double>(dual));
        }
    }
}