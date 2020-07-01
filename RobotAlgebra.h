//
// Created by sba on 26.06.20.
//

#ifndef DAK_ROBOTALGEBRA_H
#define DAK_ROBOTALGEBRA_H

#include "DualNumber.h"

// to implement:

using namespace DualNumberAlgebra;

namespace RobotAlgebra {

    template<class T>
    using Mat3 = Eigen::Matrix<T,3,3>;

    template<class T>
    using Vec3 = Eigen::Matrix<T,3,1>;

    template<class T>
    using Adj = Eigen::Matrix<T,6,6>;

    struct Distance {
        double d;
        double phi;
        Distance(double d, double phi) : d(d), phi(phi) {};
    };

    template<class T>
    void decompose(const Adj<T> & adj, Mat3<T> &rot, Vec3<T> &trans) {
        rot = adj.topLeftCorner(3,3);

        Mat3<T> sk = adj.bottomLeftCorner(3,3) * rot.inverse();
        trans << sk(2,1), sk(0,2), sk(1,0);
    }

    //generic
    template<class T>
    Adj<T> adjungate(const Mat3<T> & diag, const Mat3<T> & dual) noexcept {
        Adj<T> data;
        data.topLeftCorner(3, 3)     = diag;
        data.topRightCorner(3, 3)     = Mat3<T>::Zero(3, 3);
        data.bottomLeftCorner(3, 3)     = dual;
        data.bottomRightCorner(3, 3)     = diag;
        return data;
    }

    template<class T>
    Adj<T> waschmaschine(const Adj<T> & adj) noexcept {
        return - adj*adj;
    }

    template<class T>
    Adj<T> sandwich(const Adj<T> & adj) noexcept {
        return Adj<T>::Identity(6,6) - adj;
    }

    template<class T>
    Adj<T> adjungateDualAngle(const DualNumber<T> &dn) noexcept {
        Mat3<T> i3 = Mat3<T>::Identity(3,3);
        Mat3<T> diag = dn.getReal() * i3;
        Mat3<T> dual = dn.getDual() * i3;
        return adjungate( diag, dual);
    }

    template<class T>
    Mat3<T> cross_mat(const Vec3<T> &data) noexcept {
        Mat3<T> m;
        m << T(0), -data(2,0), data(1,0),
             data(2,0), T(0), -data(0,0),
            -data(1,0), data(0,0), T(0);
        return m;
    }

    template<class T>
    class Pluecker {
    public:
        Vec3<T> m;
        Vec3<T> n;

        Pluecker(const Vec3<T> & n, const Vec3<T> & m) noexcept {
            this->n = n;
            this->m = m;
        }

        static Pluecker<T> fromPoints(const Vec3<T> &a, const Vec3<T> &b) noexcept {
            Vec3<T> n = b - a;
            n = n / n.norm();
            Vec3<T> m = a.cross(n);
            return Pluecker<T>(n,m);
        }

        Adj<T> adjungatePluecker() const noexcept{
            return adjungate(cross_mat(n), cross_mat(m));
        }

        Vec3<DualNumber<T>> getDualLine() const noexcept{
            Vec3<DualNumber<T>> l;
            l<< DualNumber(n(0,0), m(0,0)),
                DualNumber(n(1,0), m(1,0)),
                DualNumber(n(2,0), m(2,0));
            return l;
        }

        Adj<T> REGG(DualNumber<T> phi) const noexcept {
            Adj<T> adj = this->adjungatePluecker();
            Adj<T> wm = waschmaschine(adj);
            Adj<T> sw = sandwich(wm);

            return adjungateDualAngle(cos(phi)) * wm +
                       adjungateDualAngle(sin(phi)) * adj +
                       sw;

        }

        DualNumber<T> dot(const Pluecker &rhs) const noexcept {
            return this->getDualLine().transpose() * rhs.getDualLine();
        }

        Distance getDistance(const Pluecker &rhs) const noexcept {
            DualNumber<T> dn = this->dot(rhs);
            double phi = acos(dn.getReal());
            double d = - dn.getDual() / sin(phi);
            return Distance(d, phi);
        }


    };
}

#include <iostream>
#include <iomanip>

namespace std {

    template<class T>
    std::ostream &operator<<(std::ostream &stream, RobotAlgebra::Pluecker<T> const &d) {
        return stream << d.n << std::endl << d.m << std::endl;
    }
}

#endif //DAK_ROBOTALGEBRA_H
