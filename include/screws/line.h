//
// Created by sba on 19.07.21.
//

#ifndef DUAL_ALGEBRA_KINEMATICS_LINE_H
#define DUAL_ALGEBRA_KINEMATICS_LINE_H

#include "screws/screw.h"

class UnitLine;
class DualFrame;

class Line : public Screw {
protected:
    explicit Line(const Vec6 &data) noexcept;
public:
    explicit Line(const DirectionVector &n, const PointVector &a) noexcept;
    explicit Line(const PointVector &a, const PointVector &b) noexcept;

    /**
     * \brief Align the line thus remove the pitch of the screw
     *
     * This yields to a moment vector which is orthogonal to the lines direction
     *
     * @return Aligned line
     */
    Line align() const;

    Line operator-() const noexcept;

    /**
     * \brief Normalize the screw
     *
     * @return Normalized screw
     */
    UnitLine normalize() const;

    Screw translate(double value) const noexcept;

    friend Line operator*(const DualFrame &lhs, const Line &rhs) noexcept;
};


#endif //DUAL_ALGEBRA_KINEMATICS_LINE_H
