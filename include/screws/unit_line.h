//
// Created by sba on 19.07.21.
//

#ifndef DUAL_ALGEBRA_KINEMATICS_UNIT_LINE_H
#define DUAL_ALGEBRA_KINEMATICS_UNIT_LINE_H

#include "screws/unit_screw.h"

class UnitLine : public UnitScrew {
protected:
    explicit UnitLine(const Vec6 &data) noexcept;
public:
    explicit UnitLine(const UnitDirectionVector &n, const PointVector &a) noexcept;
    explicit UnitLine(const PointVector &a, const PointVector &b) noexcept;

    /**
     * \brief Align the line thus remove the pitch of the screw
     *
     * This yields to a moment vector which is orthogonal to the lines direction
     *
     * @return Aligned line
     */
    UnitLine align() const;

    UnitLine operator-() const noexcept;

    /**
     * \brief Normalize the screw
     *
     * @return Normalized screw
     */
    UnitLine normalize() const;

    UnitScrew translate(double value) const noexcept;
    Line rotate(double value) const noexcept;
    Screw transform(DualNumberAlgebra::DualNumber value) const noexcept;

    friend UnitLine operator*(const DualFrame &lhs, const UnitLine &rhs) noexcept;
};



#endif //DUAL_ALGEBRA_KINEMATICS_UNIT_LINE_H
