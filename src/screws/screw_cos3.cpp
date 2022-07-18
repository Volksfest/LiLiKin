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

DualNumber
UnitLine::acos3(const UnitLine &a, const UnitLine &b) const noexcept {
    try {
        // The orthogonals are not computable if something is coinciding
        // In that case there is no transformation at all
        // See the catch
        UnitLine orthogonal_a = this->orthogonal(a).to_line();
        UnitLine orthogonal_b = this->orthogonal(b).to_line();

        // Check for (anti-)parallelity
        bool an_parallel = Compare::is_equal(abs(a.n() * this->n()), 1.0);
        bool bn_parallel = Compare::is_equal(abs(b.n() * this->n()), 1.0);

        // Compute the orientiation of rotation by the triple product
        // This may be 0 and thus yield to unprecise sign determination
        // But it is not critical as it also means a half-circle rotation where the sign is completely irrelevant
        double ornt = sgn( this->n() * cross(orthogonal_a.n(), orthogonal_b.n()) );
        // Compute the rotation by the acos
        auto arg = orthogonal_a.n() * orthogonal_b.n();
        // Values need to be clipped due to floating point precision
        if ( arg < -1) {arg = -1;}
        if ( arg >  1) {arg =  1;}
        // Compute the angle with its respective sign
        double angle = acos(arg ) * ornt;

        // The default translation with any parallel line is zero
        double translation = 0;
        if (!an_parallel && !bn_parallel) {
            // This computes the distance of the orthogonal containing orthogonal plane to n
            // Differently explained: This gives the constant offset of a Hesse normal form
            // The orthogonal is obviously orthogonal to the line direction and thus inside a orthogonal plane
            double plane_d_a = this->n() * orthogonal_a.get_canonical_anchor();
            double plane_d_b = this->n() * orthogonal_b.get_canonical_anchor();

            // The diffrence of the plane is the needed translation
            translation = plane_d_b - plane_d_a;
        }

        return DualNumber(
                angle,
                translation
        );
    } catch(const std::domain_error &e) {
        return DualNumber(); // 0+0ϵ
    }
}
