//
// Created by sba on 19.07.21.
//

#include <memory>

#include <util/precision.h>
#include "embedded_types/dual_frame.h"
#include "embedded_types/dual_embedded_matrix.h"
#include "screws/unit_line.h"
#include "screws/screw.h"
#include "base/vector.h"

UnitLine::UnitLine(const Vec6 &data) noexcept: Screw(data) {}

UnitLine::UnitLine(const UnitDirectionVector &n, const PointVector &a) noexcept : Screw(n, MomentVector(a.cross(n))) {}

UnitLine::UnitLine(const UnitDirectionVector &n, const MomentVector &m) : Screw(n, m) {
    if ( !Compare::is_zero(this->n() * this->m())) {
        throw std::domain_error("Moment has to be orthogonal to direction");
    }
}

UnitLine::UnitLine(const DirectionVector &n, const PointVector &a) noexcept : UnitLine(n.normal(), a) {}

UnitLine::UnitLine(const PointVector &a, const PointVector &b) : UnitLine(DirectionVector(b-a).normal(), a){}

UnitLine UnitLine::operator-() const noexcept {
    return UnitLine(-this->data);
}

UnitLine operator*(const DualFrame &lhs, const UnitLine &rhs) noexcept {
    return UnitLine(lhs.data * rhs.data);
}

Screw UnitLine::transform(DualNumberAlgebra::DualNumber value) const noexcept {
    return Screw(
            DirectionVector(this->n() * value.real()),
            MomentVector(this->n() * value.dual() + this->m() * value.real()));
}

Screw UnitLine::transform(double rotation, double translation) const noexcept {
    return this->transform(DualNumberAlgebra::DualNumber(rotation, translation));
}

std::tuple<Screw, Screw, Screw> UnitLine::decomposition(const Screw &l) const {
    auto ortho = this->orthogonal(l);
    auto rej = this->rejection(l, &ortho);
    auto proj = this->projection(l, &rej);
    return std::make_tuple(ortho,rej, proj);
}

Screw UnitLine::orthogonal(const Screw &l) const {
    Vector nm = this->n().cross(l.m()) + this->m().cross(l.n()); // na x mb + ma x nb
    try {
        return Screw(
                DirectionVector(this->n().cross(l.n())),// na x nb
                MomentVector(nm) // na x mb + ma x nb
        );
    } catch(const std::invalid_argument &) {
        try {
            // this is structural the same as above, but it needs swapping which is not possible with the types used here
            // the result could be a direction with zero-norm
            // another type with switched direction and moment would be necessary or something like that
            // but this would be much work only for this purpose.
            return Screw(
                    DirectionVector(nm), // na x mb + ma x nb
                    MomentVector(this->m().cross(l.m())) // ma x mb
                    );

        } catch(const std::invalid_argument &) {
            throw std::domain_error("Orthogonal of coinciding lines is not possible");
        }
    }
}

Screw UnitLine::rejection(const Screw &l, const Screw *orthogonal) const {
    const Screw &o = (orthogonal == nullptr) ? this->orthogonal(l) : *orthogonal;

    // "double cross product" yields to something like 180° rotation
    // (for orthogonal vectors at least. for non orthogonal it yields to a projection)
    // by negating, the rejection "goes to the same direction" as the line l except it is projected to be orthogonal to this line
    return -this->orthogonal(o);
}

Screw UnitLine::projection(const Screw &l, const Screw *rejection) const noexcept {
    try {
        // It is possible that the rejection cannot be calculated.
        // This mostly means that the lines are coinciding.
        const Screw &r = (rejection == nullptr) ? this->rejection(l) : *rejection;

        // the projection is the the part of the screw, which does not contain anyhting from the rejection.
        // regarding lines it does not make much sense but in screws the pitches are different
        auto p = l - r;
        return p;
    } catch(const std::domain_error &) {
        // It doubt the reference line itself is a good projection as without respect of the norm/pitch it is always the same line.
        return *this;
    }
}

PointVector UnitLine::point_project(const UnitLine &l) const noexcept {
    auto d = l.get_canonical_anchor() - this->get_canonical_anchor();
    // Probably, the point projection gets expensive but this makes it stable and without any consideration of cases
    auto n = this->line_cross(l);
    return PointVector(
            this->get_canonical_anchor() - (l.n().cross(l.n().cross(this->n())) * d / (n*n)) * this->n());
}

UnitLine UnitLine::parallel_through_anchor(const PointVector &new_anchor) const noexcept {
    return UnitLine(this->n(), new_anchor);
}

PointVector
UnitLine::intersect(const UnitLine &l) const {
    //use normalized screws as we don't want to have pitches and scaling by direction vector
    UnitLine a = *this;
    UnitLine b = l;

    // check if lines are coplanar and non-parallel
    if(cross(a.n(), b.n()).is_zero()) {
        throw std::domain_error("parallel lines cannot intersect");
    }
    if ( !(Compare::is_zero(a.m() * b.n() + a.n() * b.m()))) {
        throw std::domain_error("skew lines cannot intersect");
    }

    if ( a.m().is_zero()) {
        if(b.m().is_zero()) {
            // if both moments are zero, both lines are going through zero. The intersection is quite obvious, isn't it?
            return PointVector(0,0,0);
        }
        // if a's moment is zero, just swap the lines to ensure a's moment is the non-zero one
        std::swap(a,b);
    }
    auto c = cross(a.m(), b.m());
    Vector m;

    // if the moments are parallel or b's moment is zero use the moment of the orthogonal of the lines.
    // This one goes also through the intersection point is then guaranteed to be orthogonal.
    if (b.m().is_zero() || c.is_zero()) {
        // last one could be enhanced by checking if b.m is zero
        m = cross(a.m(), b.n()) + cross(a.n(), b.m());
    } else {
        m = b.m();
    }

    // Denominator can't be zero as n can't be zero and m is checked not to be zero
    // In theory taken from wikipedia but also self derived especially with origin lines
    // With lots of algebraic transformations this results from the fact, that any point on line is
    // orthogonal to the moment. Thus the intersection needs to be orthogonal to both moments and is thus the scaled cross product of both moments
    // the scaling can be derived from the algebraic substituion of the moments
    return PointVector(cross(m, a.m()) / (m * a.n()));
}

DualNumberAlgebra::DualNumber UnitLine::get_distance(const UnitLine &rhs) const noexcept {
    auto n1 = this->n();
    auto n2 = rhs.n();


    double d = 0;
    try {
        auto n12 = this->line_cross(rhs);
        d = (n1.cross(this->m()) - n2.cross(rhs.m())) * n12.normal();
    } catch(...) {
        // Do nothing
    }

    double prod = n1 * n2;
    if (prod > 1) prod = 1;
    if (prod < -1) prod = -1;

    return DualNumberAlgebra::DualNumber(
            acos(prod),
            -d
            );
}

DirectionVector UnitLine::line_cross(const UnitLine &rhs) const {
    try {
        return DirectionVector(this->n().cross(rhs.n())); // na x nb
    } catch(const std::invalid_argument &) {
        try {
            return DirectionVector(this->n().cross(rhs.m()) + this->m().cross(rhs.n())); // na x mb + ma x nb
        } catch(const std::invalid_argument &) {
            throw std::domain_error("Orthogonal direction of coinciding lines is not possible");
        }
    }
}

double UnitLine::get_distance(const PointVector &rhs) const noexcept {
    // The same as above but shortened. You can construct a parallel line through the point.
    // The rotation is zero thus it is ignored and only the dual part as a real number is given back.
    return (this->m() - cross(rhs,this->n())).norm();
}

LineRelation UnitLine::get_relation_to(const UnitLine &l) const noexcept {
    auto distance = this->get_distance(l);

    bool shared_point = false;

    if (Compare::is_zero(distance.dual())) {
        shared_point = true;
    }

    if (Compare::is_zero(distance.real())) {
        if (shared_point) {
            return LineRelation::COINCIDE;
        } else {
            return LineRelation::PARALLEL;
        }
    }

    if (Compare::is_equal(distance.real(),M_PI)) {
        if (shared_point) {
            return LineRelation::ANTI_COINCIDE;
        } else {
            return LineRelation::ANTI_PARALLEL;
        }
    }

    if (shared_point) {
        return LineRelation::INTERSECT;
    } else {
        return LineRelation::SKEW;
    }
}

UnitLine UnitLine::orthogonal_through_anchor(const PointVector &anchor) const {
    UnitLine l = this->to_line();
    try {
        DirectionVector n(
                anchor - cross(l.n(), l.m()) - l.n() * (l.n() * anchor)
        );
        return UnitLine(n, anchor);
    } catch( std::invalid_argument &) {
        throw std::domain_error("Cannot create a orthogonal through anchor if the anchor is on the line");
    }
}

PointVector UnitLine::point_project(const PointVector &p) const noexcept {
   return PointVector(this->n() * p * this->n() + this->n().cross(this->m()));
}

PointVector UnitLine::orthogonal_plane_projection(const PointVector &plane_anchor, const PointVector &point) const noexcept {
    auto diff = (this->point_project(point) - this->point_project(plane_anchor));
    double sign = (this->n() * diff) > 0 ? 1.0 : -1.0;
    return PointVector(point - diff.norm() * sign * this->n());
}

