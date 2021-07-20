//
// Created by sba on 14.07.21.
//

#ifndef DUAL_ALGEBRA_KINEMATICS_SCREW_H
#define DUAL_ALGEBRA_KINEMATICS_SCREW_H

#include <eigen3/Eigen/Eigen>

class AdjungateMatrix;

class Line;
class UnitScrew;
class DualEmbeddedMatrix;

class MomentVector;
class DirectionVector;
class PointVector;

class Projection;

namespace DualNumberAlgebra{
    class DualNumber;
}

/**
 * \brief The representation of a screw
 *
 * The screw is defined as a 6x1 Eigen Matrix.
 * It has not a normalized direction vector and has a potential pitch in its moment
 *
 * The head three elements of the data vector define the direction and the last three the moment.
 *
 */
class Screw{

protected:
    using Vec6 = Eigen::Matrix<double,6,1>;

    Vec6 data;
    explicit Screw(const Vec6 &data) noexcept;

public:
    Screw() = delete;

    /**
     * \brief Create a screw with a direction and a moment according to pluecker coordinates
     *
     * There are no conditions for the direction vector and the moment vector.
     * If the direction vector is known to be unit it should be considered to create a UnitScrew.
     * If the moment is known to be orthogonal to direction a Line should be created.
     *
     * @param n The direction vector (unit vector!)
     * @param m The moment vector according to pluecker coordinates
     */
    explicit Screw(const DirectionVector &n, const MomentVector &m) noexcept;

    /**
     * \brief Return the direction vector
     * @return The direction vector
     */
    DirectionVector n() const noexcept;

    /**
     * \brief Return the moment vector
     * @return  The moment vector
     */
    MomentVector m() const noexcept;

    /**
     * \brief Calculate the scalar dual number product between screws as an operator overload
     *
     * \todo create a free function with concrete name
     *
     *
     * @param rhs The other screw
     * @return The rotation and distance between the screws
     */
    DualNumberAlgebra::DualNumber operator*(const Screw &rhs) const noexcept;

    /**
    * \brief Calculate a norm of the skew
    *
    * This works better if the skew is interpreted as a dual vector - by having the moment vector as the dual part
    *
    * The real-part is the norm of the direction.
    * The dual-part is the scaled pitch. The scaling is given by the norm of the direction.
    *
    * This yields to a decomposition of the skew to a product by the unit line and the norm
    * such that $ = phi * Lamda_hat
    *
    * @return The norm of the skew
    */
    DualNumberAlgebra::DualNumber norm() const noexcept;

    /**
     * \brief Retrieve the canoncical anchor of the line
     *
     * In short: retrieve an arbitrarily anchor.
     *
     * @return Canonical Anchor
     */
    PointVector get_canonical_anchor() const noexcept;

    /**
     * \brief The (additive) negation and thus invert of direction of the line
     *
     * The lines are oriented as the direction vector can be negated.
     * Thus the lines are actually spears
     * \todo naming with lines and spears
     * @return The negated line
     */
    Screw operator-() const noexcept;

    /**
     * \brief The difference between two Screws
     *
     * There is probably also a complex geometry interpretation which is unimportant here
     *
     * @param rhs
     * @return
     */
    Screw operator-(const Screw &rhs) const noexcept;

    /**
     * \brief Align the line thus remove the pitch of the screw
     *
     * This yields to a moment vector which is orthogonal to the lines direction
     *
     * @return Aligned line
     */
    Line align() const;

    /**
     * \brief Normalize the screw
     *
     * @return Normalized screw
     */
    UnitScrew normalize() const;

    /**
     * \brief Create a decomposition of a line as projection to this line
     * @param l Line to decompose
     * @return The decomposition in a decomposition
     */
    Projection project(const Screw &l) const noexcept;

    /**
     * \brief Project the anchor of a line to the this line
     * @param l The line where the anchor should be projected
     * @return The projection point
     */
    PointVector point_project(const Screw &l) const noexcept;

    /**
     * \brief Create a translated line going through the new anchor
     * @param new_anchor The new anchor
     * @return The new parallel line goint through the new anchor
     */
    Line parallel_through_anchor(const PointVector &new_anchor) const noexcept;



    DualNumberAlgebra::DualNumber acos3(const Screw &a, const Screw &b) const noexcept;

    virtual bool no_rotation() const noexcept;

    friend Screw operator*(const DualEmbeddedMatrix &lhs, const Screw &rhs) noexcept;

    friend bool operator==(const Screw &lhs, const Screw &rhs);

    friend bool operator!=(const Screw &lhs, const Screw &rhs);
};

struct Projection{
    Screw p;
    Screw r;
    Screw o;
};


#endif //DUAL_ALGEBRA_KINEMATICS_SCREW_H
