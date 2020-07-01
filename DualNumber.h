//
// Created by sba on 12.06.20.
//

#ifndef DAK_DUALNUMBER_H
#define DAK_DUALNUMBER_H

#include <stdexcept>
#include <cmath>

namespace DualNumberAlgebra {

    template<class T>
    class DualNumber {
    public:
        explicit DualNumber<T>() noexcept : DualNumber<T>(0, 0) {}

        explicit DualNumber<T>(const T &real) noexcept : DualNumber<T>(real, 0) {}

        explicit DualNumber<T>(const T &real, const T &dual) noexcept {
            this->_real = real;
            this->_dual = dual;
        }


        // basic
        DualNumber<T> &operator=(const DualNumber<T> &rhs) noexcept {
            this->_real = rhs._real;
            this->_dual = rhs._dual;
            return *this;
        }

        DualNumber<T> &operator=(const T &real) noexcept {
            this->_real = real;
            this->_dual = 0;
            return *this;
        }

        DualNumber<T> operator+(const DualNumber<T> &rhs) const noexcept {
            return DualNumber<T>(this->_real + rhs._real, this->_dual + rhs._dual);
        }

        DualNumber<T> operator+(const T &real) const noexcept {
            return DualNumber<T>(this->_real + real, this->_dual);
        }

        DualNumber<T> operator+() const noexcept {
            return DualNumber<T>(this->_real, this->_dual);
        }

        DualNumber<T> operator*(const DualNumber<T> &rhs) const noexcept {
            return DualNumber<T>(this->_real * rhs._real, this->_dual * rhs._real + this->_real * rhs._dual);
        }

        DualNumber<T> operator*(const T &real) const noexcept {
            return DualNumber<T>(this->_real * real, this->_dual * real);
        }

        DualNumber<T> operator-(const DualNumber<T> &rhs) const noexcept {
            return DualNumber<T>(this->_real - rhs._real, this->_dual - rhs._dual);
        }

        DualNumber<T> operator-(const T &real) const noexcept {
            return DualNumber<T>(this->_real - real, this->_dual);
        }

        DualNumber<T> operator-() const noexcept {
            return DualNumber<T>(-this->_real, -this->_dual);
        }

        DualNumber<T> conjugate() const noexcept {
            return DualNumber<T>(this->_real, -this->_dual);
        }

        DualNumber<T> inverse() const { //
            if (this->_real == 0) {
                throw std::logic_error("Cannot invert a dual number with real part equals to zero");
            }
            return DualNumber<T>(1 / this->_real, -this->_dual / this->_real / this->_real);
        }

        T norm() const {
            return this->_real >= 0 ? this->_real : -this->_real;
        }

        DualNumber<T> operator/(const DualNumber<T> &rhs) const {
            if (this->_real == 0) {
                throw std::logic_error("Cannot divide by a dual number with real part equals to zero");
            }
            return DualNumber<T>(this->_real / rhs._real,
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

    template<class T>
    DualNumber<T> sin(const DualNumber<T> & phi) noexcept {
        const T &r = phi.getReal();
        const T &d = phi.getDual();
        return DualNumber<T>(std::sin(r), d * std::cos(r));
    }

    template<class T>
    DualNumber<T> asin(const DualNumber<T> & w) noexcept {
        const T &r = w.getReal();
        const T nr = std::asin(r);
        const T &d = w.getDual();
        return DualNumber<T>(nr, d / std::cos(nr));
    }

    template<class T>
    DualNumber<T> cos(const DualNumber<T> & phi) noexcept {
        const T &r = phi.getReal();
        const T &d = phi.getDual();
        return DualNumber<T>(std::cos(r), -d * std::sin(r));
    }

    template<class T>
    DualNumber<T> acos(const DualNumber<T> & w) noexcept {
        const T &r = w.getReal();
        const T nr = std::acos(r);
        const T &d = w.getDual();
        return DualNumber<T>(nr, -d / std::sin(nr));
    }

    template<class T>
    DualNumber<T> tan(const DualNumber<T> & phi) noexcept {
        const T &r = phi.getReal();
        const T &d = phi.getDual();
        const T cr = std::cos(r);
        return DualNumber<T>(std::tan(r), d / cr / cr);
    }

    template<class T>
    DualNumber<T> atan(const DualNumber<T> & w) noexcept {
        const T &r = w.getReal();
        const T nr = std::atan(r);
        const T &d = w.getDual();
        const T cnr = std::cos(nr);
        return DualNumber<T>(nr, -d * cnr * cnr);
    }

    template <class T>
    DualNumber<T> atan2(const DualNumber<T> &y, const DualNumber<T> &x) noexcept {
        const T &xr = x.getReal();
        const T &yr = y.getReal();
        const T &xd = x.getDual();
        const T &yd = y.getDual();
        return DualNumber<T>(std::atan2(yr, xr), (xr*yd - xd * yr) / ( xr * xr + yr * yr));
    }

    template<class T>
    inline const DualNumber<T> conj(const DualNumber<T>& x)  {
        return x.conjugate();
    }

    // for double type as lhs
    template<class T>
    DualNumber<T> operator +(const T &lhs, const DualNumber<T> &rhs) noexcept {
        return rhs + lhs;
    }

    // for constant pure dual number const defining
    // e.g 1.0 + 3.0_s can be written as a const for a dual number
    DualNumber<double> operator "" _s(long double dual) noexcept{
        return DualNumber<double>(0.0, dual);
    }

    // Simpler outputstream in the form "(a+bs)"
    #include <iostream>
    template<class T>
    std::ostream &operator<<(std::ostream &stream, DualNumberAlgebra::DualNumber<T> const &d) {
        T dual = d.getDual();
        return (stream << "(" << d.getReal() << (dual < 0 ? "-" : "+") << (dual < 0 ? -dual : dual) << "s)");
    }
}

// Eigen definition. Has to be defined for each value
namespace Eigen {
    template<>
    struct NumTraits<DualNumberAlgebra::DualNumber<double>> : NumTraits<double>
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

#endif //DAK_DUALNUMBER_H
