//
// Created by sba on 06.07.21.
//

/*
 * In the beginning, this was intended to be a template library.
 * But it is hard to write a legit library in C++ without the concepts(of C++20).
 * As this is mostly specialized for this library I dropped the idea.
 * So no templates, just pure double Dualnumbers as they are needed.
 */

#ifndef DAK_DUAL_NUMBER_H
#define DAK_DUAL_NUMBER_H

#include <stdexcept>
#include <cmath>

namespace DualNumberAlgebra {

    class DualNumber {
    private:
        double _real;
        double _dual;
        
    public:
        // basic constructors
        explicit DualNumber() noexcept : DualNumber(0, 0) {}
        explicit DualNumber(const double &real, const double &dual) noexcept : _real(real), _dual(dual) {}

        // intended to be usable as implicit constructor
        DualNumber(const double &real) noexcept : DualNumber(real, 0) {} // NOLINT

        // special constructors
        DualNumber(const DualNumber &copy) = default;
        DualNumber(DualNumber &&move) = default;

        // basic assignments
        DualNumber &operator=(const DualNumber &rhs) = default;
        DualNumber &operator=(DualNumber &&rhs) = default;
        DualNumber &operator=(const double &real) noexcept;


        ////// basic operators
        // unary
        DualNumber operator+() const noexcept;
        DualNumber operator-() const noexcept;

        // binary
        DualNumber operator+(const DualNumber &rhs) const noexcept;
        DualNumber operator*(const DualNumber &rhs) const noexcept;
        DualNumber operator-(const DualNumber &rhs) const noexcept;
        DualNumber operator/(const DualNumber &rhs) const; //does an exception: division by 0

        // compound operators
        DualNumber& operator+=(const DualNumber &rhs) noexcept;
        DualNumber& operator*=(const DualNumber &rhs) noexcept;
        DualNumber& operator-=(const DualNumber &rhs) noexcept;
        DualNumber& operator/=(const DualNumber &rhs); //does an exception: division by 0

        // more complex operations
        DualNumber conjugate() const noexcept;
        DualNumber inverse() const;
        double norm() const noexcept;
        double norm_square() const noexcept;

        bool is_zero() const noexcept;

        // projections = "getter"
        double real() const noexcept;
        double dual() const noexcept;
    };

    // mapping with derivatives
    DualNumber sin(const DualNumber & phi) noexcept;
    DualNumber asin(const DualNumber & w) noexcept;
    DualNumber cos(const DualNumber & phi) noexcept;
    DualNumber acos(const DualNumber & w) noexcept;
    DualNumber tan(const DualNumber & phi) noexcept;
    DualNumber atan(const DualNumber & w) noexcept;
    DualNumber atan2(const DualNumber &y, const DualNumber &x) noexcept;
    DualNumber sqrt(const DualNumber &x) noexcept;

    // class independant definitions for real lhs
    // binary with real number; only lhs as dual number. Rhs has to be defined afterwards
    //DualNumber operator+(DualNumber lhs, double rhs)  noexcept;
    DualNumber operator+(double  lhs, DualNumber rhs)  noexcept;
    //DualNumber operator*(DualNumber lhs, double rhs)  noexcept;
    DualNumber operator*(double  lhs, DualNumber rhs)  noexcept;
    //DualNumber operator-(DualNumber lhs, double rhs)  noexcept;
    DualNumber operator-(double  lhs, DualNumber rhs)  noexcept;
    //DualNumber operator/(DualNumber lhs, double rhs); // does an exception: division by 0
    DualNumber operator/(double  lhs, DualNumber rhs);

    bool operator==(const DualNumber a, const DualNumber b) noexcept;

    // const literals to write dual numbers more beautiful

    namespace literals {
        DualNumber operator "" _s(long double dual) noexcept;
        DualNumber operator "" _s(unsigned long long dual) noexcept;
    }

    // output stream representation
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

#endif //DAK_DUAL_NUMBER_H
