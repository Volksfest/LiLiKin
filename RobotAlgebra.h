//
// Created by sba on 26.06.20.
//

#ifndef DAK_ROBOTALGEBRA_H
#define DAK_ROBOTALGEBRA_H

#include "DualNumber.h"

using namespace DualNumberAlgebra;

namespace RobotAlgebra {

    template<class T>
    using Mat3 = Eigen::Matrix<T,3,3>;

    template<class T>
    using Vec3 = Eigen::Matrix<T,3,1>;

    template<class T>
    using Mat6 = Eigen::Matrix<T,6,6>;

    template<class T>
    using Mat4 = Eigen::Matrix<T,4,4>;

    template<class T>
    using Pluecker = Eigen::Matrix<T,6,1>;

    struct Distance {
        double d;
        double phi;
        Distance(double d, double phi) : d(d), phi(phi) {};
        Distance(const DualNumber<double> & dn) noexcept {
            phi = acos(dn.getReal());
            d = - dn.getDual() / sin(phi);
        };
    };

    template<class T>
    Mat3<T> create_cross_mat(const Vec3<T> &data) noexcept {
        Mat3<T> m;
        m << T(0), -data(2,0), data(1,0),
                data(2,0), T(0), -data(0,0),
                -data(1,0), data(0,0), T(0);
        return m;
    }

    template<class T>
    Pluecker<T> create_pluecker_from_points(const Vec3<T> &a, const Vec3<T> &b) noexcept {
        Vec3<T> n = b - a;
        n = n / n.norm();
        Vec3<T> m = a.cross(n);
        Pluecker<T> p;
        p << n, m;
        return p;
    }

    template<class T>
    Mat6<T> create_adjungate(const Mat3<T> & diag, const Mat3<T> & dual) noexcept {
        Mat6<T> data;
        data.topLeftCorner(3, 3)     = diag;
        data.topRightCorner(3, 3)     = Mat3<T>::Zero(3, 3);
        data.bottomLeftCorner(3, 3)     = dual;
        data.bottomRightCorner(3, 3)     = diag;
        return data;
    }

    template<class T>
    Mat6<T> adjungate_dualangle(const DualNumber<T> &dn) noexcept {
        Mat3<T> i3 = Mat3<T>::Identity(3,3);
        Mat3<T> diag = dn.getReal() * i3;
        Mat3<T> dual = dn.getDual() * i3;
        return create_adjungate(
                static_cast<Mat3<T>>(dn.getReal() * i3),
                static_cast<Mat3<T>>(dn.getDual() * i3));
    }

    template<class T>
    Mat6<T> adjungate_pluecker(const Pluecker<T> & p) noexcept {
        return create_adjungate(
                create_cross_mat(static_cast<Vec3<T>>(p.head(3))),
                create_cross_mat(static_cast<Vec3<T>>(p.tail(3))));
    }

    template<class T>
    Mat6<T> convert_mat4_to_mat6(const Mat4<T> & hom) noexcept {
        Vec3<T> trans = hom.topRightCorner(3,1);
        Mat3<T> R = hom.topLeftCorner(3,3);
        Mat3<T> skewR = create_cross_mat(trans) * R;
        return create_adjungate(R, skewR);
    }

    template<class T>
    Mat4<T> convert_mat6_to_mat4(const Mat6<T> & adj) noexcept {
        Mat4<T> hom = Mat4<T>::Identity(4,4);
        hom.topLeftCorner(3,3) = adj.topLeftCorner(3,3);

        Mat3<T> sk = adj.bottomLeftCorner(3,3) * hom.topLeftCorner(3,3).inverse();

        hom(0,3) = sk(2,1);
        hom(1,3) = sk(0,2);
        hom(2,3) = sk(1,0);

        return hom;
    }

    template<class T>
    Mat3<T> get_rot(const Mat4<T> & hom) noexcept {
        return hom.topLeftCorner(3,3);
    }

    template<class T>
    Vec3<T> get_trans(const Mat4<T> & hom) noexcept {
        return hom.topRightCorner(3,1);
    }

    template<class T>
    Vec3<DualNumber<T>> get_dualline(const Pluecker<T> & p) noexcept{
        Vec3<DualNumber<T>> l;
        l<< DualNumber(p(0,0), p(3,0)),
            DualNumber(p(1,0), p(4,0)),
            DualNumber(p(2,0), p(5,0));
        return l;
    }

    template<class T>
    Mat6<T> REGG(const Pluecker<T> & p, const DualNumber<T> & phi) noexcept {
        // see papers for the signs Waschmaschine, Sandwich, ...
        Mat6<T> orthoterm = adjungate_pluecker(p);
        Mat6<T> uniterm = - orthoterm * orthoterm; //eq -i^2
        Mat6<T> nullterm = Mat6<T>::Identity(6,6) - uniterm;

        return adjungate_dualangle(cos(phi)) * uniterm +
               adjungate_dualangle(sin(phi)) * orthoterm +
               nullterm;

    }

    template<class T>
    DualNumber<T> dot(const Pluecker<T> &lhs, const Pluecker<T> &rhs) noexcept {
        return lhs.get_dualLine().transpose() * rhs.get_dualLine();
    }

    // syntactic sugar for scalar product
    template<class T>
    DualNumber<T> operator*(const Pluecker<T> &lhs, const Pluecker<T> &rhs) noexcept {
        return dot(lhs, rhs);
    }

    template<class T>
    Distance getDistance(const Pluecker<T> &lhs, const Pluecker<T> &rhs) noexcept {
        DualNumber<T> dn = lhs*rhs;
        // implies that T is double! TODO type safe
        return Distance(dn);
    }
}

#endif //DAK_ROBOTALGEBRA_H
