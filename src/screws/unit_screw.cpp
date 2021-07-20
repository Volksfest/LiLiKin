//
// Created by sba on 19.07.21.
//

#include "base/vector.h"
#include "screws/unit_screw.h"
#include "screws/unit_line.h"
#include "embedded_types/dual_frame.h"

UnitScrew::UnitScrew(const Vec6 &data) noexcept : Screw(data) {}

UnitScrew::UnitScrew(const UnitDirectionVector &n, const MomentVector &m) noexcept : Screw(n,m) {}

UnitLine UnitScrew::align() const {
    return UnitLine(this->n(), this->get_canonical_anchor());
}

UnitScrew UnitScrew::operator-() const noexcept {
    return UnitScrew(-this->data);
}

UnitDirectionVector UnitScrew::n() const noexcept {
    return UnitDirectionVector(this->data(0), this->data(1), this->data(2));
}

UnitScrew UnitScrew::normalize() const {
    return Screw::normalize();
}

UnitLine UnitScrew::parallel_through_anchor(const PointVector &new_anchor) const noexcept {
    return UnitLine(this->n(), new_anchor);
}

bool UnitScrew::no_rotation() const noexcept {
    return true;
}

UnitScrew operator*(const DualFrame &lhs, const UnitScrew &rhs) noexcept {
    return UnitScrew(lhs.data * rhs.data);
}
