//
// Created by sba on 14.07.21.
//

#ifndef DUAL_ALGEBRA_KINEMATICS_SCREW_H
#define DUAL_ALGEBRA_KINEMATICS_SCREW_H

#include <eigen3/Eigen/Eigen>

class AdjungateMatrix;

class UnitLine;
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
    /**
     * \brief typedef for a Eigen 6x1 matrix, which is used as a column vector
     */
    using Vec6 = Eigen::Matrix<double,6,1>;

    /**
     * \brief Internal storage as Eigen type
     */
    Vec6 data;

    /**
     * \brief Protected constructor with raw Eigen type
     *
     * It is protected to enforce a high level user to use the other constructors.
     * Especially important to avoid zero direction vectors
     *
     * @param data Data in Eigen type
     */
    explicit Screw(const Vec6 &data) noexcept;

public:
    /**
     * \brief Deleted default constructor
     *
     * Especially to avoid zero direction vectors
     */
    Screw() = delete;

    /**
     * \brief Create a screw with a direction and a moment according to pluecker coordinates
     *
     * There are no conditions for the direction vector and the moment vector. (except n != 0)
     *
     * @param n The direction vector (anything but zero)
     * @param m The moment vector
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
     * This is just the dot product of the screw in dual vector representation.
     * The semantics are the same as the usual dual product with the cosine relation allthough it has some caveats.
     *
     * Thus use the UnitLine::get_distance(const UnitLine &) to get the corrected angular and translational distance.
     *
     * @param rhs The other screw
     * @return The dot product between screws
     */
    DualNumberAlgebra::DualNumber operator*(const Screw &rhs) const noexcept;

    /**
     * \brief Calculate a norm of the skew
     *
     * This works better if the skew is interpreted as a dual vector - by having the moment vector as the dual part
     *
     * The real-part is the norm of the direction.
     * The dual-part is the scaled moment. The scaling is given by the norm of the direction.
     *
     * This yields to a decomposition of the skew to a product by the unit line and the norm
     * such that screw = phi * unit_line
     *
     * The real-part is semantically the lenght of the screw
     * The dual-part is semantically the pitch of the screw
     *
     * There is a caveat, that the real-part is always one for UnitLine and UnitScrew allthough it should be 0.
     * This discrepancy can be detected by checking Screw::no_rotation()
     *
     * @return The norm of the skew
     */
    DualNumberAlgebra::DualNumber norm() const noexcept;

    /**
     * \brief Retrieve the canoncical anchor of the line
     *
     * Not checked, but this should be the closest point to the origin.
     *
     * Therefore canonical means here the closest to the origin which gives uniqueness.
     *
     * @return Canonical Anchor
     */
    PointVector get_canonical_anchor() const noexcept;

    /**
     * \brief The (additive) inversion of direction of the screw
     *
     * The lines are oriented as the direction vector can be negated.
     * Thus the lines are actually spears
     * @return The negated line
     */
    Screw operator-() const noexcept;

    /**
     * \brief The difference between two Screws
     *
     * There is probably also a complex geometry interpretation which is unimportant here.
     * Algebraically it is possible to compute but only if the lines are different.
     *
     * \exception std::domain_error
     * @param rhs right-hand-side
     * @return The difference screw
     */
    Screw operator-(const Screw &rhs) const;

    /**
     * \brief Align and normalize the screw.
     *
     * This yields to a moment vector which is orthogonal to the lines direction and an unit line direction.
     * This also means that there is no pitch in the screw anymore
     *
     * @return Normalized and pitchless line
     */
    UnitLine to_line() const noexcept;

    /**
     * \brief Generic Matrix-Vector product in context of embedded matrix and screws
     *
     * This is a generic Matrix-Vector without much semantics.
     *
     * @param lhs The left-hand-side matrix
     * @param rhs The right-hand-side vector
     * @return Resulting product screw
     */
    friend Screw operator*(const DualEmbeddedMatrix &lhs, const Screw &rhs) noexcept;

    /**
     * \brief Check for equality element-wise
     * @param lhs left-hand-side
     * @param rhs right-hand-side
     * @return True if equal
     */
    friend bool operator==(const Screw &lhs, const Screw &rhs);

    /**
     * \brief Check for inequality element-wise
     * @param lhs left-hand-side
     * @param rhs right-hand-side
     * @return True if unequal
     */
    friend bool operator!=(const Screw &lhs, const Screw &rhs);

    /**
     * \brief Pretty print of a screw
     * @param stream Where to write
     * @param rhs The printed screw
     * @return
     */
    friend std::ostream &operator<<(std::ostream &stream, const Screw &rhs);
};

#endif //DUAL_ALGEBRA_KINEMATICS_SCREW_H
