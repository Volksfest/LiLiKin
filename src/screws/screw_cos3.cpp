//
// Created by sba on 20.07.21.
//

#include "screws/screw.h"
#include "screws/unit_line.h"
#include "base/dual_number.h"
#include "base/vector.h"

#include "util/precision.h"

using DualNumberAlgebra::DualNumber;

double
sgn(double x) {
    return x > 0.0 ? 1.0 : -1.0;
}

// generic case
DualNumber
acos3_generic(const UnitLine &a, const UnitLine &b, const UnitLine &n) noexcept {
    auto a_pro = n.project(a);
    UnitLine rejection_a = a_pro.r.to_line();

    auto b_pro = n.project(b);
    UnitLine rejection_b = b_pro.r.to_line();

    double dot_prod = n.n() * cross(a.n(), b.n());
    double ornt = sgn(dot_prod);

    // if a and b are nearly parallel, the sign is hard to compute with the cross product
    // thus calculate the difference of the intersection which should be aligned to the normals direction
    // otherwise it is negative
    if( Compare::is_zero(dot_prod) ) {
        ornt = sgn ( (n.intersect(rejection_b) - n.intersect(rejection_a)) * n.n() );
    }
    return rejection_a.get_distance(rejection_b) * ornt;
}

//only translation or 180 Degrees
DualNumber
acos3_parallel_lines(const UnitLine &a, const UnitLine &b, const UnitLine &n) noexcept {
    auto point_projection_a = n.point_project(a);
    auto point_projection_b = n.point_project(b);

    auto projection_diff = point_projection_b - point_projection_a;

    double ornt = sgn(n.n() * projection_diff);

    // a_n and b_n are parallel, thus the product can only be +1 or -1
    // Depended on that there is a 180° rotation or not
    // the rest is just a translation as the two lines are parallel
    return DualNumber(a.n() * b.n() > 0 ? 0.0 : M_PI, projection_diff.norm() * ornt);
}

// only rotation
DualNumber
acos3_missing_rejections(const UnitLine &a, const UnitLine &b, const UnitLine &n) noexcept {
    // Just in this special case, either a or b can be coincident to n
    // in that case any kind of projection is not possible and throws an exception
    // with the exception received the coincident case is catched and no transformation (0+0ϵ) is returned
    try {
        auto a_pro_o = n.project(a).o.to_line();
        auto b_pro_o = n.project(b).o.to_line();

        double ornt = sgn(
                n.n() * cross(a_pro_o.n(), b_pro_o.n()) );

        // thanks to floating point precision the arg can be greater or lesser than 1 or -1 respectively.
        // Infortunately, the acos cannot be calculated then.
        // Thus, the argument is trimmed inside the domain
        auto arg = a_pro_o.n() * b_pro_o.n();
        if ( arg < -1) {arg = -1;}
        if ( arg >  1) {arg =  1;}
        // real acos not dual acos!
        return DualNumber(acos(arg ) * ornt, 0);
    } catch(const std::invalid_argument &e) {
        return DualNumber(); // 0+0ϵ
    }
}

DualNumber
Screw::acos3(const Screw &a, const Screw &b) const noexcept {
    UnitLine a_l = a.to_line();
    UnitLine b_l = b.to_line();
    UnitLine n_l = this->to_line();

    bool ab_parallel = Compare::is_equal(abs(a_l.n() * b_l.n()), 1.0);
    bool an_parallel = Compare::is_equal(abs(a_l.n() * n_l.n()), 1.0);
    bool bn_parallel = Compare::is_equal(abs(b_l.n() * n_l.n()), 1.0);

    // Check for "Generic case" - everthing is skewed
    if (!ab_parallel && !an_parallel && !bn_parallel) {
        return acos3_generic(a_l, b_l, n_l);

    // something is parallel
    } else {

        // lines are parallel but not to reference
        // yields in only translation
        if (!an_parallel && !bn_parallel) {
            return acos3_parallel_lines(a_l, b_l, n_l);

        // at least one line is parallel to the reference yielding to a non possible rejection
        // thus the orthogonals are used

        // disclaimer: both lines have to be parallel to the reference
        // a single line parallel to the reference simple makes no sense and will result to something unexpected
        } else {
            return acos3_missing_rejections(a_l, b_l, n_l);
        }
    }
}
