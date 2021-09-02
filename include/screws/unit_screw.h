//
// Created by sba on 19.07.21.
//

#ifndef DUAL_ALGEBRA_KINEMATICS_UNIT_SCREW_H
#define DUAL_ALGEBRA_KINEMATICS_UNIT_SCREW_H

#include "screws/screw.h"

class UnitLine;
class UnitDirectionVector;
class DualFrame;

/**
 * \brief A screw with a unit direction vector
 *
 * Semantically this describes a screw without a rotation as it only has a pitch,
 *   which can be interpreted as a translation.
 */
class UnitScrew : public Screw {
protected:
    /**
     * \brief Protected constructor with raw Eigen type
     *
     * It is protected to enforce a high level user to use the other constructor.
     * Especially important to avoid wrong direction vector format.
     *
     * @param data Data in Eigen type
     */
    explicit UnitScrew(const Vec6 &data) noexcept;
public:
    /**
     * \brief Constructor of the screw
     *
     * The direction vector has to be a unit direction vector ensuring a unit norm.
     * The moment can be arbitrary.
     * If a point is needed create a UnitLine.
     * If a point is needed but still with a non orthogonal moment, then still create a UnitLine
     *   and give it a pitch with translate.
     *
     * @param n The unit direction vector
     * @param m The moment of the screw
     */
    explicit UnitScrew(const UnitDirectionVector &n, const MomentVector &m) noexcept;

    /**
     * \brief Align the line thus remove the pitch of the screw
     *
     * This yields to a moment vector which is orthogonal to the lines direction
     *
     * @return Aligned line
     */
    UnitLine align() const;

    /**
     * \brief The (additive) inversion of direction of the unit screw
     *
     * The lines are oriented as the direction vector can be negated.
     * Thus the lines are actually spears
     * @return The negated unit screw
     */
    UnitScrew operator-() const noexcept;

    /**
     * \brief Return the direction vector as a unit direction vector
     * @return The direction vector
     */
    UnitDirectionVector n() const noexcept;

    /**
     * \brief Normalize the screw
     *
     * @return Normalized screw
     */
    UnitScrew normalize() const;

    /**
     * \brief Create a translated line going through the new anchor
     * @param new_anchor The new anchor
     * @return The new parallel line going through the new anchor
     */
    UnitLine parallel_through_anchor(const PointVector &new_anchor) const noexcept;

    /**
     * \brief Checks if the screw does not rotate
     * \todo rename
     * @return True if it has no rotation
     */
    bool no_rotation() const noexcept override;

    /**
     * \brief Apply a rotation to the unit screw be determining how much it should rotate
     *
     * This rotation makes the UnitScrew to a Screw where the norm of the direction vector determines the rotation length
     * @param value The amount of rotation in radian
     * @return The resulting screw
     */
    Screw rotate(double value) const noexcept;

    /**
     * \brief Transformation of the UnitScrew to a frame
     *
     * This should be considered as a usual coordinate transformation
     *   and thus the frame and the screw should have the same origin coordinate system.
     *
     * @param lhs The frame
     * @param rhs The screw to transform
     * @return The screw transformed by the frame
     */
    friend UnitScrew operator*(const DualFrame &lhs, const UnitScrew &rhs) noexcept;
};

#endif //DUAL_ALGEBRA_KINEMATICS_UNIT_SCREW_H
