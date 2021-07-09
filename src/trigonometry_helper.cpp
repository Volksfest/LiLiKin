//
// Created by sba on 06.07.21.
//

#include "trigonometry_helper.h"
#include <iostream>

// some partial functions. actually not necessary
double
sgn(double x) {
    return x > 0.0? 1.0 : -1.0;
}

std::vector<DualNumber>
solve_trig(const DualNumber &a,const DualNumber &b, const DualNumber &c) {

    DualNumber dd = a*a+b*b-c*c;
    DualNumber rad;

    if(dd.real() <= 0) {
        // TODO epsilon
        if ( dd.real() == 0.0) {
            return {};
        } else {
            throw std::logic_error("No solution possible");
        }
    } else {
        DualNumber d = DualNumberAlgebra::sqrt(dd);
        rad = atan2(d,c);
    }

    DualNumber pre = atan2(b,a);
    std::vector<DualNumber> solutions;

    if (rad.is_zero()) {
        solutions.push_back(pre);
    } else {
        solutions.push_back(pre + rad);
        solutions.push_back(pre - rad);
    }

    return solutions;
}

// generic case
DualNumber
acos3_generic(const Pluecker &a, const Pluecker &b, const Pluecker &n) noexcept {

    auto a_pro = n.project(a);
    Pluecker rejection_a = a_pro.r.normalize().align();

    auto b_pro = n.project(b);
    Pluecker rejection_b = b_pro.r.normalize().align();

    double ornt = sgn(
            n.n() * cross(a.n(),b.n()) );

    return acos(rejection_a * rejection_b) * ornt;

}

DualNumber
acos3_parallel_lines(const Pluecker &a, const Pluecker &b, const Pluecker &n) noexcept {
    auto point_projection_a = n.point_project(a);
    auto point_projection_b = n.point_project(b);

    auto projection_diff = point_projection_b - point_projection_a;

    double ornt = sgn(n.n() * projection_diff);

    return DualNumber(0.0, projection_diff.norm() * ornt);
}

DualNumber
acos3_missing_rejections(const Pluecker &a, const Pluecker &b, const Pluecker &n) noexcept {
    auto a_pro = n.project(a);
    auto b_pro = n.project(b);

    double ornt = sgn(
            n.n() * cross(a_pro.o.n(), b_pro.o.n()) );
    auto res = acos(a_pro.o * b_pro.o) * ornt;

    return res;
}

DualNumber
acos3(const Pluecker &a, const Pluecker &b, const Pluecker &n) noexcept {

    bool ab_parallel = abs(a.n() * b.n()) > 0.9999;
    bool an_parallel = abs(a.n() * n.n()) > 0.9999;
    bool bn_parallel = abs(b.n() * n.n()) > 0.9999;

    if (!ab_parallel && !an_parallel && !bn_parallel) {
        // generic case
        //std::cout << "Generic Case" << std::endl;
        return acos3_generic(a, b, n);

    } else {
        // something is parallel
        if (!an_parallel && !bn_parallel) {
            // lines are parallel but not to reference
            // yields in only translation
            // std::cout << "Lines are parallel" << std::endl;
            return acos3_parallel_lines(a, b, n);
        } else {
            // at least one line is parallel to the reference yielding to a non possible rejection
            // thus the orthogonals are used
            // std::cout << "Parallelity to Reference" << std::endl;
            return acos3_missing_rejections(a, b, n);
        }
    }

}