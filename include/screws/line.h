//
// Created by sba on 19.07.21.
//

#ifndef DUAL_ALGEBRA_KINEMATICS_LINE_H
#define DUAL_ALGEBRA_KINEMATICS_LINE_H

#include "screws/screw.h"

class UnitLine;
class DualFrame;

/**
 * \brief A screw without a pitch
 *
 * This can also be interpreted as a screw, which only describes a rotation around the direction vector around its norm.
 */
class Line : public Screw {
protected:
    /**
     * \brief Protected constructor with raw Eigen type
     *
     * It is protected to enforce a high level user to use the other constructor.
     * Especially important to avoid wrong direction vector format.
     *
     * @param data Data in Eigen type
     */
    explicit Line(const Vec6 &data) noexcept;
public:
    /**
     * \brief Construction of the line with direction and anchor point
     *
     * @param n The direction vector
     * @param a The anchor point
     */
    explicit Line(const DirectionVector &n, const PointVector &a) noexcept;

    /**
     * \brief Construction of the line with two points
     *
     * The direction goes from the frist to the second point
     *
     * \exception std::logic_error Thrown by DirectionVector if a and b are equal
     * @param a The first point
     * @param b The second point
     */
    explicit Line(const PointVector &a, const PointVector &b);

    /**
     * \brief Align the line thus remove the pitch of the screw
     *
     * This yields to a moment vector which is orthogonal to the lines direction.
     *
     * In this case it is overwritten to do nothing.
     *
     * @return Aligned line
     */
    Line align() const;

    /**
     * \brief The (additive) inversion of direction of the line
     *
     * The lines are oriented as the direction vector can be negated.
     * Thus the lines are actually spears
     * @return The negated line
     */
    Line operator-() const noexcept;

    /**
     * \brief Normalize the line
     *
     * @return Normalized unit line
     */
    UnitLine normalize() const;

    /**
     * \brief Apply a translation to the line
     *
     * This translation makes the Line to a Screw where the norm of the direction vector determines the rotation length
     *  and the value the scaled pitch
     * @param value The amount of translation. The unit is given as "length unit" and depends on how the user defines his point vectors.
     * @return The resulting screw
     */
    Screw translate(double value) const noexcept;

    /**
     * \brief Transformation of the Line to a frame
     *
     * This should be considered as a usual coordinate transformation
     *   and thus the frame and the line should have the same origin coordinate system.
     *
     * @param lhs The frame
     * @param rhs The line to transform
     * @return The line transformed by the frame
     */
    friend Line operator*(const DualFrame &lhs, const Line &rhs) noexcept;
};


#endif //DUAL_ALGEBRA_KINEMATICS_LINE_H
