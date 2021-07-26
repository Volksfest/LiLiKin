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

Screw::Screw(const Vec6 &data) noexcept : data(data) {}

Screw::Screw(const DirectionVector &n, const MomentVector &m) noexcept : data(Vec6::Zero(6,1)) {
    this->data << n.data, m.data;
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

Screw Screw::operator-(const Screw &rhs) const noexcept {
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

    std::unique_ptr<Screw> op;
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
    auto p = l - r;

    return {p,r,o};
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
    return lhs.data == rhs.data && lhs.no_rotation() == rhs.no_rotation();
}

bool operator!=(const Screw &lhs, const Screw &rhs) {
    return ! (lhs == rhs);
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
    if (!std::isfinite(angle.dual())) {
        return DualNumberAlgebra::DualNumber(
                // still copy the real part as it can be either 0 or Pi (depending on the direction of the spear)
                angle.real(),
                // as the lines are parallel the canonical anchors should be on the orthogonal plane through zero. their distance is also the distance between the lines.
                (t.get_canonical_anchor() - r.get_canonical_anchor()).norm());
    } else {
        return angle;
    }
}
