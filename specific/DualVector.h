//
// Created by sba on 29.06.20.
//

#ifndef DAK_DUALVECTOR_H
#define DAK_DUALVECTOR_H

#include "DualNumber.h"
#include "MatDefinitions.h"

namespace DualNumberAlgebra {

    class DualVector {
    private:
        Vec3 _real;
        Vec3 _dual;

    public:
        explicit DualVector() = default;

        explicit DualVector(const Vec3 &real) noexcept;

        explicit DualVector(const Vec3 &real, const Vec3 &dual) noexcept;

        // basic
        DualVector &operator=(const DualVector &rhs) = default;

        DualVector &operator=(const Vec3 &real) noexcept;

        DualVector operator+(const DualVector &rhs) const noexcept;

        DualVector operator+(const Vec3 &real) const noexcept;

        DualVector operator+() const noexcept;

        DualVector operator-(const DualVector &rhs) const noexcept;

        DualVector operator-(const Vec3 &real) const noexcept;

        DualVector operator-() const noexcept;

        DualVector conjugate() const noexcept;

        double norm() const noexcept;

        // projections = "getter"
        Vec3 getReal() const noexcept;

        Vec3 getDual() const noexcept;

        // Dot product

        DualNumber operator*(const DualVector &rhs) const noexcept;

        DualNumber operator*(const Vec3 &real) const noexcept;

    };
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
#endif //DAK_DUALVECTOR_H
