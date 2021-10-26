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
 * \brief Gives the parallelity
 */
enum Parallelity{
    SKEW, ///< Skewed lines
    INTERSECT, ///< Intersecting lines
    PARALLEL, ///< Parallel lines
    ANTI_PARALLEL, ///< Parallel lines with opposite orientation
    COINCIDE, ///< Coinciding lines
    ANTI_COINCIDE ///< Coinciding lines with opposite orientation
};

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
     * Thus use the Screw::get_distance(const Screw &) to get the corrected angular and translational distance.
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
     * \brief The (additive) inversion of direction of the line
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
     * \brief Create a decomposition of a screw as projections to this screw as reference
     * \todo determine the exception
     * @param l Screw to decompose
     * @return The decomposition
     */
    Projection project(const Screw &l) const;

    /**
     * \brief Project a screw to this screw as reference
     *
     * This is actually just the point on this screw which is closest to the other screw.
     *
     * @param l Screw which should be projected
     * @return The projection point
     */
    PointVector point_project(const Screw &l) const noexcept;

    /**
     * \brief Orthogonal projection of a point to the screw
     * @param p Point to project
     * @return The projected point
     */
    PointVector point_project(const PointVector &p) const noexcept;

    /**
     * \brief Create a translated line going through the new anchor
     * @param new_anchor The new anchor
     * @return The new parallel line going through the new anchor
     */
    Screw parallel_through_anchor(const PointVector &new_anchor) const noexcept;

    /**
     * \brief Create a line being orthogonal and going through the anchor
     * To determine a unique orthogonal line the anchor must not be on this screw
     * \exception std::invalid_argument Cannot determine an orthogonal as the anchor is on the line
     * @param anchor The anchor where the orthogonal line should go through
     * @return The orthogonal line through the anchor
     */
    Screw orthogonal_through_anchor(const PointVector &anchor) const;

    /**
     * \brief Project a point orthogonally to a plane spanned by the direction of this screw and an anchor point
     *
     * The direction of this screw given by Screw::n() is used as the hesse normal form of the plane.
     * The anchor point determines the position of the plane.
     *
     * @param plane_anchor Point where the plane goes through
     * @param point Point to project
     * @return The projected point
     */
    PointVector orthogonal_plane_projection(const PointVector &plane_anchor, const PointVector &point) const noexcept;

    /**
     * \brief Get the intersection of two screws
     * \exception std::domain_error If the screws don't intersect. Coinciding screws also throws an exception
     * @param l right-hand-side
     * @return The intersection point
     */
    PointVector intersect(const Screw &l) const;

    /**
     * \brief Get the distance between the point and the screw
     * @param rhs The evaluated point
     * @return The distance to the evaluated point
     */
    double get_distance(const PointVector &rhs) const noexcept;

    /**
     * \brief Get the distance between two screws
     * This is a corrected value. Four classes can be distinguished for a dual number \f$ x = r + \epsilon d \f$
     *
     * r = d = 0 : Coinciding screws
     *
     * r = 0 , d != 0: Parallel screws
     *
     * r != 0 , d = 0: Intersecting screws
     *
     * r != 0, d != 0: skewed screws
     *
     * Caveat: In theory, r = pi also results in coinciding and parallel screws.
     * The screws don't have the same orientation anymore but are still parallel or coinciding without orientation.
     * This can yield to many problems, so it has to be checked. See also Screw::is_parallel() and Parallelity.
     *
     * @param rhs right-hand-side
     * @return Angular and translational distance as dual number
     */
    DualNumberAlgebra::DualNumber get_distance(const Screw &rhs) const noexcept;

    /**
     * \brief Retrieve the relationship between two screws
     * \todo Rename the function
     * @see Parallelity
     * @param l
     * @return
     */
    Parallelity is_parallel(const Screw &l) const noexcept;

    /**
     * \brief Calculate the transformation of this screw to get the projection of a onto b
     *
     * This can be seen as the solver for b = exp(phi * l) a
     *
     * Where phi should be found and l is this line.
     *
     * @param a The screw a
     * @param b The screw b
     * @return The found solution for the transformation around this line
     */
    DualNumberAlgebra::DualNumber acos3(const Screw &a, const Screw &b) const noexcept;

    /**
     * \brief Checks if the screw does not rotate
     * \todo rename
     * @return True if it has no rotation
     */
    virtual bool no_rotation() const noexcept;

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
};

/**
 * \brief Struct which contains the projections
 */
struct Projection{
    /**
     * \brief Projected part of the screw
     *
     * This screw coincides with the reference screw
     */
    Screw p;

    /**
     * \brief Rejected part of the screw
     *
     * This screw is orthogonal to the reference but not to the original screw.
     * The sum of the rejection and the projection yields to the original screw again.
     */
    Screw r;

    /**
     * \brief Orthogonal part of the screw
     *
     * This screw is orthogonal to the original screw and the reference.
     * Thus is can be used to transform the original screw to the reference.
     */
    Screw o;
};


#endif //DUAL_ALGEBRA_KINEMATICS_SCREW_H
