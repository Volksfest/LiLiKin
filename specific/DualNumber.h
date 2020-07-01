//
// Created by sba on 29.06.20.
//

#ifndef DAK_DUALNUMBER_H
#define DAK_DUALNUMBER_H

#include <stdexcept>
#include <cmath>

namespace DualNumberAlgebra {

    class DualNumber {
    private:
        double _real;
        double _dual;
        
    public:
        explicit DualNumber() noexcept : DualNumber(0, 0) {}

        explicit DualNumber(const double &real) noexcept : DualNumber(real, 0) {}

        explicit DualNumber(const double &real, const double &dual) noexcept;

        // basic
        DualNumber &operator=(const DualNumber &rhs) = default;

        DualNumber &operator=(const double &real) noexcept;

        DualNumber operator+(const DualNumber &rhs) const noexcept;

        DualNumber operator+(const double &real) const noexcept;

        DualNumber operator+() const noexcept;

        DualNumber operator*(const DualNumber &rhs) const noexcept;

        DualNumber operator*(const double &real) const noexcept;

        DualNumber operator-(const DualNumber &rhs) const noexcept;

        DualNumber operator-(const double &real) const noexcept;

        DualNumber operator-() const noexcept;

        DualNumber conjugate() const noexcept;

        DualNumber inverse() const;

        double norm() const noexcept;

        DualNumber operator/(const DualNumber &rhs) const;

        // projections = "getter"
        double getReal() const noexcept;

        double getDual() const noexcept;

        // mapping with derivatices

        DualNumber sin(const DualNumber & phi) const noexcept;

        DualNumber asin(const DualNumber & w) const noexcept;

        DualNumber cos(const DualNumber & phi) const noexcept;

        DualNumber acos(const DualNumber & w) const noexcept;

        DualNumber tan(const DualNumber & phi) const noexcept;

        DualNumber atan(const DualNumber & w) const noexcept;

        DualNumber atan2(const DualNumber &y, const DualNumber &x) const noexcept;
    };

    // class independant definitions

    DualNumber operator +(const double &lhs, const DualNumber &rhs) noexcept;

    DualNumber operator "" _s(long double dual) noexcept;

    /*
    inline const DualNumber<double> conj(const DualNumber<double>& x)  { return x.conjugate(); }
    inline const DualNumber<double> real(const DualNumber<double>& x)  { return DualNumber<double>(x.getReal(),0); }
    inline DualNumber<double> imag(const DualNumber<double> &x)    { return DualNumber<double>(x.getDual(),0); }
    inline DualNumber<double> abs(const DualNumber<double>&  x)  { return DualNumber<double>(x.norm(),0); }
    inline DualNumber<double> abs2(const DualNumber<double>& x)  { return DualNumber<double>(x.norm()*x.norm(),0); }
    */

}

// Simpler outputstream in the form "(a+bs)"
#include <iostream>
namespace std {
    std::ostream &operator<<(std::ostream &stream, DualNumberAlgebra::DualNumber const &d);
}

#include <eigen3/Eigen/Eigen>

// Eigen definition. Has to be defined for each value
namespace Eigen {
    template<>
    struct NumTraits<DualNumberAlgebra::DualNumber> : NumTraits<double>
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
