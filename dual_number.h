//
// Created by sba on 12.06.20.
//

#ifndef DAK_DUAL_NUMBER_H
#define DAK_DUALNUMBER_H

#include <stdexcept>
#include <cmath>
#include <iostream>
#include <eigen3/Eigen/Eigen>

namespace DualNumberAlgebra {

    template<class T>
    class DualNumber_template {
    public:
        explicit DualNumber_template<T>() noexcept : DualNumber_template<T>(0, 0) {}

        explicit DualNumber_template<T>(const T &real) noexcept : DualNumber_template<T>(real, 0) {}

        explicit DualNumber_template<T>(const T &real, const T &dual) noexcept {
            this->_real = real;
            this->_dual = dual;
        }

        bool isZero() const noexcept {
            return this->_real == 0 && this->_dual == 0;
        }

        // basic
        DualNumber_template<T> &operator=(const DualNumber_template<T> &rhs) noexcept {
            this->_real = rhs._real;
            this->_dual = rhs._dual;
            return *this;
        }

        DualNumber_template<T> &operator=(const T &real) noexcept {
            this->_real = real;
            this->_dual = 0;
            return *this;
        }

        DualNumber_template<T> operator+(const DualNumber_template<T> &rhs) const noexcept {
            return DualNumber_template<T>(this->_real + rhs._real, this->_dual + rhs._dual);
        }

        DualNumber_template<T> operator+(const T &real) const noexcept {
            return DualNumber_template<T>(this->_real + real, this->_dual);
        }

        DualNumber_template<T> operator+() const noexcept {
            return DualNumber_template<T>(this->_real, this->_dual);
        }

        DualNumber_template<T> operator*(const DualNumber_template<T> &rhs) const noexcept {
            return DualNumber_template<T>(this->_real * rhs._real, this->_dual * rhs._real + this->_real * rhs._dual);
        }

        DualNumber_template<T> operator*(const T &real) const noexcept {
            return DualNumber_template<T>(this->_real * real, this->_dual * real);
        }

        DualNumber_template<T> operator-(const DualNumber_template<T> &rhs) const noexcept {
            return DualNumber_template<T>(this->_real - rhs._real, this->_dual - rhs._dual);
        }

        DualNumber_template<T> operator-(const T &real) const noexcept {
            return DualNumber_template<T>(this->_real - real, this->_dual);
        }

        DualNumber_template<T> operator-() const noexcept {
            return DualNumber_template<T>(-this->_real, -this->_dual);
        }

        DualNumber_template<T> conjugate() const noexcept {
            return DualNumber_template<T>(this->_real, -this->_dual);
        }

        DualNumber_template<T> inverse() const { //
            if (this->_real == 0) {
                throw std::logic_error("Cannot invert a dual number with real part equals to zero");
            }
            T inv = 1 / this->_real;
            return DualNumber_template<T>(inv, -this->_dual * inv * inv);
        }

        double norm() const {
            return static_cast<double>(this->_real >= 0 ? this->_real : -this->_real);
        }

        DualNumber_template<T> operator/(const DualNumber_template<T> &rhs) const {
            if (rhs._real == 0) {
                throw std::logic_error("Cannot divide by a dual number with real part equals to zero");
            }
            return DualNumber_template<T>(this->_real / rhs._real,
                                 this->_dual / rhs._real - this->_real * rhs._dual / rhs._real / rhs._real);
        }

        // projections = "getter"
        T getReal() const noexcept {
            return this->_real;
        }

        T getDual() const noexcept {
            return this->_dual;
        }

    private:
        T _real;
        T _dual;
    };

    // To get the conjugate operator as a static function for eigen
    template<class T>
    inline const DualNumber_template<T> conj(const DualNumber_template<T>& x)  {
        return x.conjugate();
    }

    // for double type as lhs
    template<class T>
    DualNumber_template<T> operator +(const T &lhs, const DualNumber_template<T> &rhs) noexcept {
        return rhs + lhs;
    }

    // for constant pure dual number const defining
    // e.g 1.0 + 3.0_s can be written as a const for a dual number
    DualNumber_template<double> operator "" _s(long double dual) noexcept;

    // Convenience for real part for dual angles
    double operator "" _d(long double degree) noexcept;

    // Simpler outputstream in the form "(a+bs)"

    template<class T>
    std::ostream &operator<<(std::ostream &stream, DualNumberAlgebra::DualNumber_template<T> const &d) {
        T dual = d.getDual();
        return (stream << "(" << d.getReal() << (dual < 0 ? "-" : "+") << (dual < 0 ? -dual : dual) << "s)");
    }

    DualNumberAlgebra::DualNumber_template<double> sqrt(const DualNumberAlgebra::DualNumber_template<double> & phi) noexcept;

    DualNumberAlgebra::DualNumber_template<double> sin(const DualNumberAlgebra::DualNumber_template<double> & phi) noexcept;

    DualNumberAlgebra::DualNumber_template<double> asin(const DualNumberAlgebra::DualNumber_template<double> & w) noexcept;

    DualNumberAlgebra::DualNumber_template<double> cos(const DualNumberAlgebra::DualNumber_template<double> & phi) noexcept;

    DualNumberAlgebra::DualNumber_template<double> acos(const DualNumberAlgebra::DualNumber_template<double> & w) noexcept;

    DualNumberAlgebra::DualNumber_template<double> tan(const DualNumberAlgebra::DualNumber_template<double> & phi) noexcept;

    DualNumberAlgebra::DualNumber_template<double> atan(const DualNumberAlgebra::DualNumber_template<double> & w) noexcept;

    DualNumberAlgebra::DualNumber_template<double> atan2(const DualNumberAlgebra::DualNumber_template<double> &y, const DualNumberAlgebra::DualNumber_template<double> &x) noexcept;
}

// Eigen definition. Has to be defined for each value
namespace Eigen {
    template<>
    struct NumTraits<DualNumberAlgebra::DualNumber_template<double>> : NumTraits<double>
    {
        typedef double Real;
        typedef double NonInteger;
        typedef double Nested;

        static inline Real epsilon() { return 1e-8; }

        enum {
            IsComplex = 1,
            IsInteger = 0,
            IsSigned = 1,
            RequireInitialization = 1,
            ReadCost = 1,
            AddCost = 2,
            MulCost = 4
        };
    };
}

#endif //DAK_DUAL_NUMBER_H
