//
// Created by sba on 19.07.21.
//

#include <memory>
#include <cmath>

#include "embedded_types/dual_embedded_matrix.h"

#include "screws/screw.h"
#include "screws/unit_screw.h"
#include "screws/unit_line.h"
#include "screws/line.h"

#include "base/vector.h"
#include "base/matrix3.h"
#include "base/dual_number.h"

#include "util/precision.h"

Screw::Screw(const Vec6 &data) noexcept : data(data) {}

Screw::Screw(const DirectionVector &n, const MomentVector &m) noexcept : data(Vec6::Zero(6,1)) {
    this->data << n.get(), m.get();
}

DirectionVector Screw::n() const noexcept {
    return DirectionVector(this->data(0), this->data(1), this->data(2));
}

MomentVector Screw::m() const noexcept {
    return MomentVector(this->data(3), this->data(4), this->data(5));
}

DualNumberAlgebra::DualNumber Screw::operator*(const Screw &rhs) const noexcept {
    using DualNumberAlgebra::DualNumber;
    Eigen::Matrix<DualNumber,3,1> l, r;

    l<< DualNumber(this->data(0,0), this->data(3,0)),
            DualNumber(this->data(1,0), this->data(4,0)),
            DualNumber(this->data(2,0), this->data(5,0));

    r<< DualNumber(rhs.data(0,0), rhs.data(3,0)),
            DualNumber(rhs.data(1,0), rhs.data(4,0)),
            DualNumber(rhs.data(2,0), rhs.data(5,0));

    return l.transpose() * r;
}

DualNumberAlgebra::DualNumber Screw::norm() const noexcept {
    auto r = this->n().norm();
    auto d = this->n() * this->m() / r;

    return DualNumberAlgebra::DualNumber(r,d);
}

PointVector Screw::get_canonical_anchor() const noexcept {
    auto norm = this->n().norm();
    return PointVector(this->n().cross(this->m()) / (norm * norm));
}

Screw Screw::operator-() const noexcept {
    return Screw(-this->data);
}

Screw Screw::operator-(const Screw &rhs) const {
    if (this->n() == rhs.n()) {
        throw std::domain_error("Resulting direction would be zero");
    }
    return Screw(this->data - rhs.data);
}

Line Screw::align() const {
    return Line(this->n(), this->get_canonical_anchor()); // TODO could be made faster by directly computing
}

UnitScrew Screw::normalize() const {
    return UnitScrew(this->n().normal(), MomentVector(this->m()/ this->n().norm()));
}

Projection Screw::project(const Screw &l) const {
    DualEmbeddedMatrix adj(
            SkewMatrix(this->n()),
            SkewMatrix(this->m()));

    std::unique_ptr<Screw> op(nullptr);
    if (abs(this->n() * l.n()) > 0.9999) {
        auto anchor_n = this->get_canonical_anchor();
        auto c = anchor_n - l.get_canonical_anchor();

        auto n = cross(cross(this->n(), c), this->n());

        op = std::make_unique<Screw>(Line(DirectionVector(n),
                 PointVector(anchor_n)));
    } else {
        op = std::make_unique<Screw>(adj * l);
    }
    auto o = *op;
    auto r = -adj * o;
    try {
        auto p = l - r;
        return {p,r,o};
    } catch(const std::domain_error &) {
        return {*this,r,o};
    }
}

PointVector Screw::point_project(const Screw &l) const noexcept {
    PointVector point_projection = this->get_canonical_anchor();
    double cross_norm = cross(l.n(), this->n()).norm();
    if (cross_norm != 0) {
        auto projection = SkewMatrix(l.n().normal());
        point_projection += this->n() * (projection * projection * this->n() * (this->get_canonical_anchor() - l.get_canonical_anchor()) /
                                         cross_norm / cross_norm);
    }
    return point_projection;
}

Line Screw::parallel_through_anchor(const PointVector &new_anchor) const noexcept {
    return Line(this->n(), new_anchor);
}

bool Screw::no_rotation() const noexcept {
    return false;
}

Screw operator*(const DualEmbeddedMatrix &lhs, const Screw &rhs) noexcept {
    return Screw(lhs.data * rhs.data);
}

bool operator==(const Screw &lhs, const Screw &rhs) {
    return lhs.data.isApprox(rhs.data, Compare::instance().get_precision()) && lhs.no_rotation() == rhs.no_rotation();
}

bool operator!=(const Screw &lhs, const Screw &rhs) {
    return ! (lhs == rhs);
}

PointVector
Screw::intersect(const Screw &l) const {
    //use normalized screws as we don't want to have pitches and scaling by direction vector
    UnitLine a = this->align().normalize();
    UnitLine b = l.align().normalize();

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

DualNumberAlgebra::DualNumber Screw::get_distance(const Screw &rhs) const noexcept {

    // The screws needs to be aligned and normalized otherwise the dual inner product gets scaled
    UnitLine t = this->align().normalize();
    UnitLine r = rhs.align().normalize();

    auto prod = t * r;

    // fix some precission error which results in domain errors
    // the norm of the dot products of unit vectors cannot be greater than 1
    if (prod.real() > 1.0) {
        prod = DualNumberAlgebra::DualNumber(1.0, prod.dual());
    }
    if (prod.real() < -1.0) {
        prod = DualNumberAlgebra::DualNumber(-1.0, prod.dual());
    }
    auto angle = acos(prod);

    // additional check as the translation of parallel lines cannot be expressed as the dual of a dual inner product
    // this is a degradition due to the sinus of phi being zero
    if (Compare::is_equal(abs(prod.real()), 1.0)) {
        return DualNumberAlgebra::DualNumber(
                // still copy the real part as it can be either 0 or Pi (depending on the direction of the spear)
                angle.real(),
                // as the lines are parallel the moments give the distance
                prod.real() > 0 ? (t.m() - r.m()).norm() : (t.m() + r.m()).norm() );
    } else {
        return angle;
    }
}

double Screw::get_distance(const PointVector &rhs) const noexcept {
    UnitLine t = this->align().normalize();

    // The same as above but shortened. You can construct a parallel line through the point.
    // The rotation is zero thus it is ignored and only the dual part as a real number is given back.
    return (t.m() - cross(rhs,t.n())).norm();
}

Parallelity Screw::is_parallel(const Screw &l) const noexcept {
    auto distance = this->get_distance(l);

    bool shared_point = false;

    if (Compare::is_zero(distance.dual())) {
        shared_point = true;
    }

    if (Compare::is_zero(distance.real())) {
        if (shared_point) {
            return Parallelity::COINCIDE;
        } else {
            return Parallelity::PARALLEL;
        }
    }

    if (Compare::is_equal(distance.real(),M_PI)) {
        if (shared_point) {
            return Parallelity::ANTI_COINCIDE;
        } else {
            return Parallelity::ANTI_PARALLEL;
        }
    }

    if (shared_point) {
        return Parallelity::INTERSECT;
    } else {
        return Parallelity::SKEW;
    }
}

Line Screw::orthogonal_through_anchor(const PointVector &anchor) const {
    UnitLine l = this->align().normalize();
    try {
        DirectionVector n(
            anchor - cross(l.n(), l.m()) - l.n() * (l.n() * anchor)
        );
        return Line(n, anchor);
    } catch( std::domain_error &e) {
        throw std::domain_error("Cannot create a orthogonal through anchor if the anchor is on the line");
    }
}

PointVector Screw::point_project(const PointVector &p) const noexcept {
    UnitLine l = this->align().normalize();

    return PointVector(this->n() * p * this->n() + this->n().cross(this->m()));
}

PointVector Screw::orthogonal_plane_projection(const PointVector &plane_anchor, const PointVector &point) const noexcept {
    auto diff = (this->point_project(point) - this->point_project(plane_anchor));
    double sign = (this->n() * diff) > 0 ? 1.0 : -1.0;
    return PointVector(point - diff.norm() * sign * this->n());
}
