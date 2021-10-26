//
// Created by sba on 19.07.21.
//

#include <util/precision.h>
#include "embedded_types/dual_frame.h"
#include "embedded_types/dual_embedded_matrix.h"
#include "screws/unit_line.h"
#include "screws/screw.h"
#include "base/vector.h"

UnitLine::UnitLine(const Vec6 &data) noexcept: Screw(data) {}

UnitLine::UnitLine(const UnitDirectionVector &n, const PointVector &a) noexcept : Screw(n, MomentVector(a.cross(n))) {}

UnitLine::UnitLine(const UnitDirectionVector &n, const MomentVector &m) noexcept : Screw(n, m) {
    if ( !Compare::is_zero(this->n() * this->m())) {
        throw std::domain_error("Moment has to be orthogonal to direction");
    }
}

UnitLine::UnitLine(const DirectionVector &n, const PointVector &a) noexcept : UnitLine(n.normal(), a) {}

UnitLine::UnitLine(const PointVector &a, const PointVector &b) : UnitLine(DirectionVector(b-a).normal(), a){}

UnitLine UnitLine::align() const {
    return *this;
}

UnitLine UnitLine::operator-() const noexcept {
    return UnitLine(-this->data);
}

UnitLine UnitLine::normalize() const {
    return *this;
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

bool UnitLine::no_rotation() const noexcept {
    return true;
}