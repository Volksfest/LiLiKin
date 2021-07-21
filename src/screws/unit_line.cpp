//
// Created by sba on 19.07.21.
//

#include "embedded_types/dual_frame.h"
#include "embedded_types/dual_embedded_matrix.h"
#include "screws/unit_line.h"
#include "screws/unit_screw.h"
#include "screws/line.h"
#include "screws/screw.h"
#include "base/vector.h"

UnitLine::UnitLine(const Vec6 &data) noexcept: UnitScrew(data) {}

UnitLine::UnitLine(const UnitDirectionVector &n, const PointVector &a) noexcept : UnitScrew(n, MomentVector(a.cross(n))) {}

UnitLine::UnitLine(const PointVector &a, const PointVector &b) noexcept : UnitLine(DirectionVector(b-a).normal(), a){}

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

UnitScrew UnitLine::translate(double value) const noexcept {
    return UnitScrew(this->n(), MomentVector(this->m() + value * this->n()));
}

Line UnitLine::rotate(double value) const noexcept {
    return Line(DirectionVector(value * this->n()), this->get_canonical_anchor()); // TODO optimize?
}

Screw UnitLine::transform(DualNumberAlgebra::DualNumber value) const noexcept {
    return Screw(
            DirectionVector(this->n() * value.real()),
            MomentVector(this->n() * value.dual() + this->m() * value.real()));
}
