//
// Created by sba on 19.07.21.
//

#include <cmath>

#include "embedded_types/dual_embedded_matrix.h"

#include "screws/screw.h"
#include "screws/unit_line.h"

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

UnitLine Screw::to_line() const noexcept {
    return UnitLine(this->n().normal(), this->get_canonical_anchor());
}

Screw operator*(const DualEmbeddedMatrix &lhs, const Screw &rhs) noexcept {
    return Screw(lhs.data * rhs.data);
}

bool operator==(const Screw &lhs, const Screw &rhs) {
    return lhs.data.isApprox(rhs.data, Compare::instance().get_precision());
}

bool operator!=(const Screw &lhs, const Screw &rhs) {
    return ! (lhs == rhs);
}
