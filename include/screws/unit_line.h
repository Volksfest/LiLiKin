//
// Created by sba on 19.07.21.
//

#ifndef DUAL_ALGEBRA_KINEMATICS_UNIT_LINE_H
#define DUAL_ALGEBRA_KINEMATICS_UNIT_LINE_H

#include "screws/unit_screw.h"

/**
 * \brief A screw without a pitch and a unit direction
 *
 * This can also be interpreted as a screw, without any transformation.
 * So this is a pure pluecker line as it should be defined.
 *
 */
class UnitLine : public UnitScrew {
protected:
    /**
     * \brief Protected constructor with raw Eigen type
     *
     * It is protected to enforce a high level user to use the other constructor.
     * Especially important to avoid wrong direction vector format.
     *
     * @param data Data in Eigen type
     */
    explicit UnitLine(const Vec6 &data) noexcept;
public:
    /**
    * \brief Construction of the unit line with unit direction and anchor point
    *
    * @param n The unit direction vector
    * @param a The anchor point
    */
    explicit UnitLine(const UnitDirectionVector &n, const PointVector &a) noexcept;

    /**
     * \brief Construction of the line with two points
     *
     * The direction goes from the frist to the second point.
     * The direction will be automatically normalized
     *
     * \exception std::logic_error Thrown by UnitDirectionVector if a and b are equal
     * @param a The first point
     * @param b The second point
     */
    explicit UnitLine(const PointVector &a, const PointVector &b);

    /**
     * \brief Align the line thus remove the pitch of the screw
     *
     * This yields to a moment vector which is orthogonal to the lines direction
     * In this case it is overwritten to do nothing
     *
     * @return Aligned line
     */
    UnitLine align() const;

    /**
     * \brief The (additive) inversion of direction of the line
     *
     * The lines are oriented as the direction vector can be negated.
     * Thus the lines are actually spears
     * @return The negated line
     */
    UnitLine operator-() const noexcept;

    /**
     * \brief Normalize the line
     *
     * In this case it is overwritten to do nothing
     *
     * @return Normalized screw
     */
    UnitLine normalize() const;

    /**
     * \brief Apply a translation to the line
     *
     * This translation makes the Line to a UnitScrew giving it a pitch.
     * @param value The pitch size
     * @return The resulting unit screw
     */
    UnitScrew translate(double value) const noexcept;

    /**
     * \brief Apply a rotation to the unit line be determining how much it should rotate
     *
     * This rotation makes the UnitLine to a Line where the norm of the direction vector determines the rotation length
     * @param value The amount of rotation in radian
     * @return The resulting Line
     */
    Line rotate(double value) const noexcept;

    /**
     * \brief Apply a rotation and a transloation to the unit line
     *
     * The transformation is given by a dual number, where the real part corresponds to the rotation and the dual
     *   part to the translation.
     * This yields to a screw, which actually represents a transformation
     * @param value The value of the transformation with the real-part for rotation and the dual-part for translation
     * @return The resulting Screw
     */
    Screw transform(DualNumberAlgebra::DualNumber value) const noexcept;

    /**
     * \brief Transformation of the UnitLine to a frame
     *
     * This should be considered as a usual coordinate transformation
     *   and thus the frame and the screw should have the same origin coordinate system.
     *
     * @param lhs The frame
     * @param rhs The line to transform
     * @return The line transformed by the frame
     */
    friend UnitLine operator*(const DualFrame &lhs, const UnitLine &rhs) noexcept;
};



#endif //DUAL_ALGEBRA_KINEMATICS_UNIT_LINE_H
