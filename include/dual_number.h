//
// Created by sba on 06.07.21.
//

/*
 * In the beginning, this was intended to be a template library.
 * But it is hard to write a legit library in C++ without the concepts(of C++20).
 * As this is mostly specialized for this library I dropped the idea.
 * So no templates, just pure double Dual numbers as they are needed.
 */

#ifndef DAK_DUAL_NUMBER_H
#define DAK_DUAL_NUMBER_H

#include <stdexcept>
#include <cmath>

/**
 * \brief Everything regarding the dual numbers
 *
 * Contains a lot of symbols for the algebraic structure of dual numbers
 */
namespace DualNumberAlgebra {

    /**
     * \brief An implementation for dual numbers
     *
     */
    class DualNumber {
    private:
        double _real; //!< Real part of the dual number
        double _dual; //!< Dual part of the dual number
        
    public:
        /**
         * \brief Default constructor for zero
         */
        explicit DualNumber() noexcept : DualNumber(0, 0) {}

        /**
         * \brief Default constructor to create a dual number
         * @param real Real part
         * @param dual Dual part
         */
        explicit DualNumber(const double &real, const double &dual) noexcept : _real(real), _dual(dual) {}

        /**
         * \brief Implicit conversion constructor for a double as dual number with zero dual part
         * @param real
         */
        DualNumber(const double &real) noexcept : DualNumber(real, 0) {} // NOLINT

        /**
         * \brief Default copy constructor
         * @param copy
         */
        DualNumber(const DualNumber &copy) = default;

        /**
         * \brief Default move constructor
         * @param move
         */
        DualNumber(DualNumber &&move) = default;

        /**
         * \brief Default copy assign operator
         * @param rhs
         * @return
         */
        DualNumber &operator=(const DualNumber &rhs) = default;

        /**
         * \brief Default move assign operator
         * @param rhs
         * @return
         */
        DualNumber &operator=(DualNumber &&rhs) = default;

        /**
         * \brief Conversion operator for implicit usage of doubles with zero dual part
         * @param real
         * @return
         */
        DualNumber &operator=(const double &real) noexcept;

        /**
         * \brief Just for completeness
         * @return
         */
        DualNumber operator+() const noexcept;

        /**
         * \brief Additive inverse of the dual number
         * @return
         */
        DualNumber operator-() const noexcept;

        /**
         * \brief Add two dual numbers
         * @param rhs
         * @return
         */
        DualNumber operator+(const DualNumber &rhs) const noexcept;

        /**
         * \brief Multiply two dual numbers
         * @param rhs
         * @return
         */
        DualNumber operator*(const DualNumber &rhs) const noexcept;

        /**
         * \brief Subtract two dual numbers
         * @param rhs
         * @return
         */
        DualNumber operator-(const DualNumber &rhs) const noexcept;

        /**
         * \brief Divide two dual numbers
         * \except std::locig_error If divided by zero real part
         * @param rhs
         * @return
         */
        DualNumber operator/(const DualNumber &rhs) const; //does an exception: division by 0

        /**
         * \brief Add a dual number to this
         * @param rhs
         * @return
         */
        DualNumber& operator+=(const DualNumber &rhs) noexcept;

        /**
         * \brief Multiply a dual number to this
         * @param rhs
         * @return
         */
        DualNumber& operator*=(const DualNumber &rhs) noexcept;

        /**
         * \brief Subtract a dual number to this
         * @param rhs
         * @return
         */
        DualNumber& operator-=(const DualNumber &rhs) noexcept;

        /**
         * \brief Divide a dual number to this
         * \except std::logic_error If divided by real part zero
         * @param rhs
         * @return
         */
        DualNumber& operator/=(const DualNumber &rhs); //does an exception: division by 0

        /**
         * \brief Conjugate dual number
         * So negate the dual part
         * @return
         */
        DualNumber conjugate() const noexcept;

        /**
         * \brief Invert dual number
         * \except std::logic_error If real part is zero
         * @return
         */
        DualNumber inverse() const;

        /**
         * \brief The norm of the dual number
         * @return
         */
        double norm() const noexcept;

        /**
         * \brief The square of the norm
         *
         * This is cheaper than square the norm as the norm is created by the root of this
         * @return
         */
        double norm_square() const noexcept;

        /**
         * \brief Check if dual number is zero
         * @return
         */
        bool is_zero() const noexcept;

        /**
         * \brief Project the real part
         *
         * Actually just retrieve the real part
         * @return
         */
        double real() const noexcept;

        /**
         * \brief Project the dual part
         *
         * Actually just retrieve the dual part
         * @return
         */
        double dual() const noexcept;
    };

    /**
     * \brief The sin of a dual number
     * @param phi
     * @return
     */
    DualNumber sin(const DualNumber & phi) noexcept;
    /**
     * \brief The asin of a dual number
     * @param w
     * @return
     */
    DualNumber asin(const DualNumber & w) noexcept;
    /**
     * \brief The cos of a dual number
     * @param phi
     * @return
     */
    DualNumber cos(const DualNumber & phi) noexcept;
    /**
     * \brief The acos of a dual number
     * @param w
     * @return
     */
    DualNumber acos(const DualNumber & w) noexcept;
    /**
     * \brief The tan of a dual number
     * @param phi
     * @return
     */
    DualNumber tan(const DualNumber & phi) noexcept;
    /**
     * \brief The atan of a dual number
     * @param w
     * @return
     */
    DualNumber atan(const DualNumber & w) noexcept;
    /**
     * \brief The atan2 of a dual number
     * @param y
     * @param x
     * @return
     */
    DualNumber atan2(const DualNumber &y, const DualNumber &x) noexcept;
    /**
     * \brief The square root of a dual number
     * @param x
     * @return
     */
    DualNumber sqrt(const DualNumber &x) noexcept;

    /**
     * \brief Commutative symmetry for double
     * * \see DualNumber::operator+
     * @param lhs
     * @param rhs
     * @return
     */
    DualNumber operator+(double  lhs, DualNumber rhs)  noexcept;
    /**
     * \brief Commutative symmetry for double
     * * \see DualNumber::operator*
     * @param lhs
     * @param rhs
     * @return
     */
    DualNumber operator*(double  lhs, DualNumber rhs)  noexcept;
    /**
     * \brief Commutative symmetry for double
     * \see DualNumber::operator-
     * @param lhs
     * @param rhs
     * @return
     */
    DualNumber operator-(double  lhs, DualNumber rhs)  noexcept;
    /**
     * \brief Commutative symmetry for double
     * \except std::logic_error If real part of right-hand-side is zero
     * \see DualNumber::operator/
     * @param lhs
     * @param rhs
     * @return
     */
    DualNumber operator/(double  lhs, DualNumber rhs);

    /**
     * \brief Simple comparison
     *
     * a and b has to be equal
     * @param a
     * @param b
     * @return
     */
    bool operator==(const DualNumber &a, const DualNumber &b) noexcept;

    // const literals to write dual numbers more beautiful

    /**
     * \brief Literals for simple dual number const with dual parts
     *
     * Like (3 + 5_s) for DualNumber(3,5)
     */
    namespace literals {
        DualNumber operator "" _s(long double dual) noexcept;
        DualNumber operator "" _s(unsigned long long dual) noexcept;
    }

    /**
     * \brief Implementation to display the dual number in the standard output
     * @param stream
     * @param d
     * @return
     */
    std::ostream &operator<<(std::ostream &stream, DualNumberAlgebra::DualNumber const &d);
}

#include <eigen3/Eigen/Eigen>

// Eigen definition. Has to be defined for each value
namespace Eigen {
    /**
     * \brief Eigen3 Description for the DualNumber
     */
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
