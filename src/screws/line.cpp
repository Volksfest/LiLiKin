//
// Created by sba on 19.07.21.
//

#include "base/vector.h"
#include "screws/line.h"
#include "screws/unit_line.h"
#include "embedded_types/dual_frame.h"

Line::Line(const Vec6 &data) noexcept: Screw(data) {}

Line::Line(const DirectionVector &n, const PointVector &a) noexcept : Screw(n, MomentVector(a.cross(n))) {}

Line::Line(const PointVector &a, const PointVector &b) noexcept : Screw(DirectionVector(b-a), MomentVector(a.cross(b))) {}

Line Line::align() const {
    return *this;
}

Line Line::operator-() const noexcept {
    return Line(-this->data);
}

UnitLine Line::normalize() const {
    return UnitLine(this->n().normal(), this->get_canonical_anchor()); // TODO optimize
}

Line operator*(const DualFrame &lhs, const Line &rhs) noexcept {
    return Line(lhs.data * rhs.data);
}