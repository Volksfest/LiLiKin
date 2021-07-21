//
// Created by sba on 20.07.21.
//

#include "screws/screw.h"
#include "screws/line.h"
#include "screws/unit_line.h"
#include "base/dual_number.h"
#include "base/vector.h"

using DualNumberAlgebra::DualNumber;

double
sgn(double x) {
    return x > 0.0 || abs(x) < 0.0000001 ? 1.0 : -1.0;
}

// generic case
DualNumber
acos3_generic(const UnitLine &a, const UnitLine &b, const UnitLine &n) noexcept {

    // TODO what about using orthogonals instead of rejections?

    auto a_pro = n.project(a);
    UnitLine rejection_a = a_pro.r.normalize().align();

    auto b_pro = n.project(b);
    UnitLine rejection_b = b_pro.r.normalize().align();

    if (abs(rejection_a.n() * rejection_b.n()) > 0.999999999) {
        return DualNumber(
                sgn(rejection_a.n() * rejection_b.n()) > 0 ? 0 : M_PI,
                acos(rejection_a * b_pro.o.normalize().align()).dual()
                );
    } else {
        double ornt = sgn(
                n.n() * cross(a.n(), b.n()));

        return acos(rejection_a * rejection_b) * ornt;
    }
}

//only translation
DualNumber
acos3_parallel_lines(const UnitLine &a, const UnitLine &b, const UnitLine &n) noexcept {
    auto point_projection_a = n.point_project(a);
    auto point_projection_b = n.point_project(b);

    auto projection_diff = point_projection_b - point_projection_a;

    double ornt = sgn(n.n() * projection_diff);

    return DualNumber(0.0, projection_diff.norm() * ornt);
}

// only rotation
DualNumber
acos3_missing_rejections(const UnitLine &a, const UnitLine &b, const UnitLine &n) noexcept {
    auto a_pro_o = n.project(a).o.normalize().align();
    auto b_pro_o = n.project(b).o.normalize().align();

    double ornt = sgn(
            n.n() * cross(a_pro_o.n(), b_pro_o.n()) );
    // real acos not dual acos!
    return DualNumber(acos(a_pro_o.n() * b_pro_o.n() ) * ornt, 0);
}

DualNumber
Screw::acos3(const Screw &a, const Screw &b) const noexcept {

    bool ab_parallel = abs(a.n() * b.n()) > 0.9999;
    bool an_parallel = abs(a.n() * this->n()) > 0.9999;
    bool bn_parallel = abs(b.n() * this->n()) > 0.9999;

    UnitLine a_l = a.normalize().align();
    UnitLine b_l = b.normalize().align();
    UnitLine n_l = this->normalize().align();

    if (!ab_parallel && !an_parallel && !bn_parallel) {
        // generic case
        //std::cout << "Generic Case" << std::endl;
        return acos3_generic(a_l, b_l, n_l);

    } else {
        // something is parallel
        if (!an_parallel && !bn_parallel) {
            // lines are parallel but not to reference
            // yields in only translation
            // std::cout << "Lines are parallel" << std::endl;
            return acos3_parallel_lines(a_l, b_l, n_l);
        } else {
            // at least one line is parallel to the reference yielding to a non possible rejection
            // thus the orthogonals are used
            // std::cout << "Parallelity to Reference" << std::endl;
            return acos3_missing_rejections(a_l, b_l, n_l);
        }
    }

}