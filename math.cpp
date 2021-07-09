//
// Created by sba on 21.09.20.
//

#include "RobotAlgebraTypes.h"
#include "math.h"


#include <iostream>

DualVector toDualLine(const Pluecker &p) noexcept {
    DualVector l;
    l<< DualNumber(p.data(0,0), p.data(3,0)),
            DualNumber(p.data(1,0), p.data(4,0)),
            DualNumber(p.data(2,0), p.data(5,0));
    return l;
}

double sgn(double a) {
    if (a > 0) {
        return 1;
    }
    return -1;
    /*
    if (a < 0) {
        return -1;
    }
    return 0;*/
}

DualNumber dot(const Pluecker &lhs, const Pluecker &rhs) noexcept {
    return toDualLine(lhs).transpose() * toDualLine(rhs);
    //return DualNumber_t(lhs.v() * rhs.v(), lhs.v() * rhs.w() + lhs.w() * rhs.v());
}

int solve_trig(const DualNumber &a, const DualNumber &b, const DualNumber &c, DualNumber &phi_1, DualNumber &phi_2) noexcept {


    DualNumber dd = a*a+b*b-c*c;

    DualNumber rad;

    if(dd.getReal() <= 0) {
        // TODO epsilon
        if ( dd.getReal() == 0.0) {
            return -1;
        } else {
            return 0;
        }
    } else {
        DualNumber d = sqrt(dd);
        rad = atan2(d,c);
    }

    DualNumber pre = atan2(b,a);

    if (rad.isZero()) {
        phi_1 = pre;
        return 1;
    }

    phi_1 = pre + rad;
    phi_2 = pre - rad;
    return 2;

    /*
    if (rad.getReal(2) == 0 && rad.getDual() == 0) {
        phi_1 = atan2(-b, c + a) * 2;
        return 1;
    }

    phi_1 = atan2(-b + d, c + a) * 2;
    phi_2 = atan2(-b - d, c + a) * 2;
    return 2;
     */
}

DualNumber acos3(const Pluecker &a, const Pluecker &b, const Pluecker &n) noexcept {

    /*
    Eigen::Matrix<double,3,3> mat; mat
        << a.v().data, b.v().data, n.v().data;
    */

    bool ab_parallel = abs(a.v() * b.v()) > 0.9999;
    bool an_parallel = abs(a.v() * n.v()) > 0.9999;
    bool bn_parallel = abs(b.v() * n.v()) > 0.9999;

    if (!ab_parallel && !an_parallel && !bn_parallel) {
        // generic case
        std::cout << "Generic Case" << std::endl;

        Pluecker p,r,o;

        a.project(n, p, r, o);
        Pluecker rejection_a = r.normalize().align();

        b.project(n, p, r, o);
        Pluecker rejection_b = r.normalize().align();

        double ornt = sgn(a.v().cross(b.v()) * n.v());

        return acos(dot(rejection_a, rejection_b)) * ornt;

    } else {
        // something is parallel

        if ( !an_parallel && !bn_parallel) {
            // lines are parallel but not to reference
            // yields in only translation

            std::cout << "Lines are parallel" << std::endl;

            auto anchor_a = a.v().cross(a.w());
            auto anchor_b = b.v().cross(b.w());
            auto anchor_n = n.v().cross(n.w());

            auto point_projection_a = anchor_n;
            auto point_projection_b = anchor_n;

            double cross_norm = a.v().cross(n.v()).norm();
            if (cross_norm != 0) {
                auto projection = RobotAlgebra::SkewMatrix::from(a.v());
                point_projection_a = point_projection_a +
                                     n.v() * (projection * projection * n.v() * (anchor_a - anchor_n) / cross_norm /
                                              cross_norm);
            }

            cross_norm = b.v().cross(n.v()).norm();
            if (cross_norm != 0) {
                auto projection = RobotAlgebra::SkewMatrix::from(b.v());
                point_projection_b = point_projection_b +
                                     n.v() * (projection * projection * n.v() * (anchor_b - anchor_n) / cross_norm /
                                              cross_norm);
            }

            auto projection_diff = point_projection_b - point_projection_a;
            double ornt = sgn(projection_diff * n.v());

            return DualNumber(0.0, ornt * projection_diff.norm());
        } else {
            // at least one line is parallel to the reference yielding to a non possible projection
            // thus the orthogonals are used

            std::cout << "Parallelity to Reference" << std::endl;

            //ignore p and r
            Pluecker p,r,orthogonal_a, orthogonal_b;

            a.project(n, p, r, orthogonal_a);
            b.project(n, p, r, orthogonal_b);

            double ornt = sgn(orthogonal_a.v().cross(orthogonal_b.v()) * n.v());
            auto res = acos(dot(orthogonal_a, orthogonal_b)) * ornt;

            return res;
        }
    }
}

bool solve(
        const AdjungateMatrix &ik,
        const AdjungateMatrix &zp,
        const Pluecker &l12_r,
        const Pluecker &l23_r,
        const Pluecker &l34_r,
        DualNumber &phi_1,
        DualNumber &phi_2,
        DualNumber &phi_3, bool other_solution) noexcept {



    Pluecker l12 = l12_r.normalize().align();
    Pluecker l23 = l23_r.normalize().align();
    Pluecker l34 = l34_r.normalize().align();

    RobotAlgebra::GenericMatrix6 s(ik.data * zp.data.inverse());

    RobotAlgebra::GenericMatrix6 crossterm = RobotAlgebra::GenericMatrix6::from(l23);
    RobotAlgebra::GenericMatrix6 uniterm = - crossterm * crossterm;
    RobotAlgebra::GenericMatrix6 squareterm(Eigen::Matrix<double,6,6>::Identity(6,6) - uniterm.data);

    DualNumber a( dot(l12, uniterm * l34));
    DualNumber b( dot(l12, crossterm * l34));
    DualNumber c( dot(l12, (s - squareterm) * l34));

    DualNumber phi_2a, phi_2b;
    std::cout << "Solutions of: " << a << " " << b << " " << c << std::endl;
    int sol = solve_trig(a,b,c, phi_2a, phi_2b);
    std::cout << "Solutions ("<< sol <<"): " << phi_2a << " " << phi_2b << std::endl;
    if ( sol == 0) {
        return false;
    }

    if (sol == -1) {
        auto p = Pluecker(s.data * l12.data);
        auto d = l12.v().cross(l12.w() - p.w()).norm();

        if (other_solution) {
            phi_2 = DualNumber(0, -d);
        } else {
            phi_2 = DualNumber(0, d);
        }
    } else {
        if (other_solution && sol == 2) {
            phi_2 = phi_2b;
        } else {
            phi_2 = phi_2a;
        }
    }

    auto m2 = l23.transform(phi_2);

    phi_1 = acos3(Pluecker(s.data * l34.data), Pluecker(m2.data * l34.data), Pluecker(-l12.data));
    phi_3 = acos3(Pluecker(s.data.inverse() * l12.data), Pluecker(m2.data.inverse() * l12.data), l34);

    // TODO change hack
    if (sol == -1) {
        auto intermediate_line = Pluecker(l12.transform(phi_1).data * m2.data * l23.data);
        auto offset = acos3(
                intermediate_line,
                Pluecker::fromDirection(intermediate_line.v(), HomMatrix::from(AdjungateMatrix(s.data)).p()),
                l34);
        std::cout << offset << std::endl;
        phi_3 = phi_3 - offset;
    }

    return true;

}