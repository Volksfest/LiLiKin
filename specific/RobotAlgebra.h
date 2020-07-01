//
// Created by sba on 29.06.20.
//

#ifndef DAK_ROBOTALGEBRA_H
#define DAK_ROBOTALGEBRA_H

#include "DualNumber.h"
#include "DualVector.h"

#include "MatDefinitions.h"

namespace RobotAlgebra {

    class Vector3;

    class Matrix3 {
    public:
        Mat3 data;

        Matrix3() noexcept;

        Matrix3(const Vec3 &a, const Vec3 &b, const Vec3 &c) noexcept;

        Matrix3(const Vector3 &a, const Vector3 &b, const Vector3 &c) noexcept;

        explicit Matrix3(const Mat3 &mat) noexcept;
    };

    class Vector3 {
    public:
        Vec3 data;

        Vector3() noexcept;

        Vector3(const double &a, const double &b, const double &c) noexcept;


        explicit Vector3(const Vec3 &a) noexcept;


        Matrix3 cross_mat() const noexcept;
    };

    class Adjungate {
    private:
        Adj data;
    public:
        Adjungate () = default;

        Adjungate(const Vector3 &n, const Vector3 &m) noexcept;
    };

    class Pluecker {
    public:
        Vector3 m;
        Vector3 n;

        Pluecker(const Vector3 & n, const Vector3 & m) noexcept;

        static Pluecker fromPoints(const Vector3 &a, const Vector3 &b);
/*
        DualNumberAlgebra::DualNumber<Vector3<T>> getDualLine() {
            return DualNumberAlgebra::DualNumber<Vector3<T>>(n, m);
        }*/
    };
}

#include <iostream>

namespace std {

    std::ostream &operator<<(std::ostream &stream, RobotAlgebra::Matrix3 const &d);

    std::ostream &operator<<(std::ostream &stream, RobotAlgebra::Vector3 const &d);

    std::ostream &operator<<(std::ostream &stream, RobotAlgebra::Pluecker const &d);
/*
    template<>
    std::ostream &operator<<(std::ostream &stream, DualNumberAlgebra::DualNumber<AdjungateAlgebra::Vector3<double>> const &v) {
        auto r = v.getReal();
        auto d = v.getDual();
        return stream << std::setfill(' ') << std::setw(8) << r.data(0,0) << "    " << d.data(0,0) << std::endl <<
                      std::setfill(' ') << std::setw(8) << r.data(1,0) << " +s " << d.data(1,0) << std::endl <<
                      std::setfill(' ') << std::setw(8) << r.data(2,0) << "    " << d.data(2,0) << std::endl;
    }*/
}

#endif //DAK_ROBOTALGEBRA_H
