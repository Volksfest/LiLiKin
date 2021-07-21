//
// Created by sba on 19.07.21.
//

#ifndef DUAL_ALGEBRA_KINEMATICS_UNIT_SCREW_H
#define DUAL_ALGEBRA_KINEMATICS_UNIT_SCREW_H

#include "screws/screw.h"

class UnitLine;
class UnitDirectionVector;
class DualFrame;

class UnitScrew : public Screw {
protected:
    explicit UnitScrew(const Vec6 &data) noexcept;
public:
    explicit UnitScrew(const UnitDirectionVector &n, const MomentVector &m) noexcept;

    /**
     * \brief Align the line thus remove the pitch of the screw
     *
     * This yields to a moment vector which is orthogonal to the lines direction
     *
     * @return Aligned line
     */
    UnitLine align() const;

    UnitScrew operator-() const noexcept;

    UnitDirectionVector n() const noexcept;

    /**
     * \brief Normalize the screw
     *
     * @return Normalized screw
     */
    UnitScrew normalize() const;

    UnitLine parallel_through_anchor(const PointVector &new_anchor) const noexcept;

    virtual bool no_rotation() const noexcept override;

    Screw rotate(double value) const noexcept;

    friend UnitScrew operator*(const DualFrame &lhs, const UnitScrew &rhs) noexcept;
};

#endif //DUAL_ALGEBRA_KINEMATICS_UNIT_SCREW_H
