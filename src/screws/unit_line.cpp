//
// Created by sba on 19.07.21.
//

#include "embedded_types/dual_frame.h"
#include "screws/unit_line.h"
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
