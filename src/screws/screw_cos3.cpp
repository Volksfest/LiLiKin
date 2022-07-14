//
// Created by sba on 20.07.21.
//

#include "screw.h"
#include "unit_line.h"
#include "dual_number.h"
#include "vector.h"

#include "precision.h"

using DualNumberAlgebra::DualNumber;

double
sgn(double x) {
    return x > 0.0 ? 1.0 : -1.0;
}

// generic case
DualNumber
acos3_generic(const UnitLine &a, const UnitLine &b, const UnitLine &n) noexcept {
    // As nothing is parallel to n the rejections definitively exist
    UnitLine rejection_a = n.rejection(a).to_line();
    UnitLine rejection_b = n.rejection(b).to_line();

    // Compute the orientiation of rotation by the triple product
    // This may be 0 and thus yield to unprecise sign determination
    // But it is not critical as it also means a half-circle rotation where the direction is completely irrelevant
    double ornt_angle = sgn( n.n() * cross(a.n(), b.n()) );
    // This computes the distance of the rejection containing orthogonal plane to n
    double plane_d_a = n.n() * rejection_a.get_canonical_anchor();
    double plane_d_b = n.n() * rejection_b.get_canonical_anchor();
    // The sign of the difference gives us the direction of the translation
    double ornt_trans = sgn( plane_d_b - plane_d_a);

    // Compute the dual angle
    // Unfortunately, the sign of the dual part is relatively unusable
    auto dual_angle = rejection_a.get_distance(rejection_b);

    // Thus, create a new dual angle where the individual signs utilized
    // The dual part always is signed and thus the absolute value is used
    return DualNumber(
            ornt_angle * dual_angle.real(),
            ornt_trans * abs(dual_angle.dual())
            );
}

// only rotation
DualNumber
acos3_missing_rejections(const UnitLine &a, const UnitLine &b, const UnitLine &n) noexcept {
    // If either a or b is coincident to n the orthogonal is also non-computable
    // In that case any kind of meaningful projection is not possible and throws an exception
    // With the exception received the coincident case is catched and no transformation (0+0ϵ) is returned (as the lines are coinciding already)
    try {
        auto a_pro_o = n.orthogonal(a).to_line();
        auto b_pro_o = n.orthogonal(b).to_line();

        double ornt = sgn(
                n.n() * cross(a_pro_o.n(), b_pro_o.n()) );

        // thanks to floating point precision the arg can be greater or lesser than 1 or -1 respectively.
        // Unfortunately, the acos cannot be calculated then.
        // Thus, the argument is trimmed inside the domain
        auto arg = a_pro_o.n() * b_pro_o.n();
        if ( arg < -1) {arg = -1;}
        if ( arg >  1) {arg =  1;}
        // real acos not dual acos!
        return DualNumber(acos(arg ) * ornt, 0);
    } catch(const std::domain_error &e) {
        return DualNumber(); // 0+0ϵ
    }
}

DualNumber
UnitLine::acos3(const UnitLine &a, const UnitLine &b) const noexcept {
    bool an_parallel = Compare::is_equal(abs(a.n() * this->n()), 1.0);
    bool bn_parallel = Compare::is_equal(abs(b.n() * this->n()), 1.0);

    if (!an_parallel && !bn_parallel) {
        return acos3_generic(a, b, *this);

        // at least one line is parallel to the reference yielding to a non possible rejection
        // thus the orthogonals are used instead
    } else {
        return acos3_missing_rejections(a, b, *this);
    }
}
