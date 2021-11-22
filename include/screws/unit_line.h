//
// Created by sba on 19.07.21.
//

#ifndef DUAL_ALGEBRA_KINEMATICS_UNIT_LINE_H
#define DUAL_ALGEBRA_KINEMATICS_UNIT_LINE_H

#include "screws/screw.h"

/**
 * \brief Description of the relationship between lines
 */
enum LineRelation {
    SKEW, ///< Skewed lines
    INTERSECT, ///< Intersecting lines
    PARALLEL, ///< Parallel lines
    ANTI_PARALLEL, ///< Parallel lines with opposite orientation
    COINCIDE, ///< Coinciding lines
    ANTI_COINCIDE ///< Coinciding lines with opposite orientation
};

class UnitDirectionVector;

class DualFrame;

/**
 * \brief A screw without a pitch and a unit direction
 *
 * This can also be interpreted as a screw, without any transformation.
 * So this is a pure pluecker line as it should be defined.
 *
 */
class UnitLine : public Screw {
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
    * \brief Construction of the unit line with unit direction and moment
    * \exception std::domain_error Moment is not orthogonal to direction
    *
    * @param n The unit direction vector
    * @param a The anchor point
    */
    explicit UnitLine(const UnitDirectionVector &n, const MomentVector &m);

    /**
    * \brief Construction of the unit line with direction and anchor point
    *
    * The direction will be normalized
    *
    * @param n The direction vector
    * @param a The anchor point
    */
    explicit UnitLine(const DirectionVector &n, const PointVector &a) noexcept;

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
     * \brief The (additive) inversion of direction of the line
     *
     * The lines are oriented as the direction vector can be negated.
     * Thus the lines are actually spears
     * @return The negated line
     */
    UnitLine operator-() const noexcept;

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
     * \brief Apply a rotation and a translation to the unit line
     *
     * The transformation is given by a dual number, where the real part corresponds to the rotation and the dual
     *   part to the translation.
     * This yields to a screw, which actually represents a transformation
     * @param value The value of the transformation with the real-part for rotation and the dual-part for translation
     * @return The resulting Screw
     */
    Screw transform(double rotation, double translation) const noexcept;

    /**
     * \brief Create the orthogonal between two screws
     *
     * \exception std::domain_error If the lines are coinciding, there is no orthogonal.
     * @param l Screw to find the common orthogonal to
     * @return The common orthogonal
     */
    Screw orthogonal(const Screw &l) const;

    /**
     * \brief Create the rejection of a screw in respect to this line
     *
     * \exception std::domain_error If the lines are coinciding, there is no rejection.
     * @param l Screw to find the rejection with respect to this line
     * @return The rejection of l
     */
    Screw rejection(const Screw &l, const Screw *orthogonal = nullptr) const;

    /**
     * \brief Create the projection of a screw in respect to this line
     *
     * @param l Screw to find the projection
     * @return The projection of l
     */
    Screw projection(const Screw &l, const Screw *rejection = nullptr) const noexcept;

    /**
     * \brief Create a decomposition of a line as projections to this line as reference
     *
     * This is a helper method for retrieving all projections needed for the decomposition.
     * \exception std::domain_error If the lines are coinciding, there is no orthogonal.
     * \see UnitLine::orthogonal UnitLine::rejection UnitLine::projection
     * @param l Line to decompose
     * @return The decomposition in the order Orthogonal, Rejection, Projection
     */
    std::tuple<Screw, Screw, Screw> decomposition(const Screw &l) const;

    /**
     * \brief Project a line to this line as reference
     *
     * This is actually just the point on this line which is closest to the other line.
     * This also means, if the projections to each other are the same, the lines intersects in that point
     *
     * @param l Screw which should be projected
     * @return The projection point
     */
    PointVector point_project(const UnitLine &l) const noexcept;

    /**
     * \brief Orthogonal projection of a point to the line
     * @param p Point to project
     * @return The projected point
     */
    PointVector point_project(const PointVector &p) const noexcept;

    /**
     * \brief Create a translated line going through the new anchor
     * @param new_anchor The new anchor
     * @return The new parallel line going through the new anchor
     */
    UnitLine parallel_through_anchor(const PointVector &new_anchor) const noexcept;

    /**
     * \brief Create a line being orthogonal and going through the anchor
     * To determine a unique orthogonal line the anchor must not be on this screw
     * \exception std::invalid_argument Cannot determine an orthogonal as the anchor is on the line
     * @param anchor The anchor where the orthogonal line should go through
     * @return The orthogonal line through the anchor
     */
    UnitLine orthogonal_through_anchor(const PointVector &anchor) const;

    /**
     * \brief Project a point orthogonally to a plane spanned by the direction of this line and an anchor point
     *
     * The direction of this line given by UnitLine::n() is used as the normal direction of the plane.
     * The anchor point determines the position of the plane.
     *
     * @param plane_anchor Point where the plane goes through
     * @param point Point to project
     * @return The projected point
     */
    PointVector orthogonal_plane_projection(const PointVector &plane_anchor, const PointVector &point) const noexcept;

    /**
     * \brief Get the intersection of two lines
     * \exception std::domain_error If the lines don't intersect. Coinciding lines also throw an exception as they don't have a single intersection point.
     * @param l right-hand-side
     * @return The intersection point
     */
    PointVector intersect(const UnitLine &l) const;

    /**
     * \brief Get the distance between the point and the line
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
     * This can yield to many problems, so it has to be checked. See also UnitLine::get_relation_to() and LineRelation.
     *
     * @param rhs right-hand-side
     * @return Angular and translational distance as dual number
     */
    DualNumberAlgebra::DualNumber get_distance(const UnitLine &rhs) const noexcept;

    /**
     *
     */
    DirectionVector line_cross(const UnitLine &rhs) const;

    /**
     * \brief Retrieve the relationship between two screws
     * @see LineRelation
     * @param l
     * @return
     */
    LineRelation get_relation_to(const UnitLine &l) const noexcept;

    /**
     * \brief Calculate the transformation of this line to get the projection of a onto the projection of b
     *
     *
     * @param a The line a
     * @param b The line b
     * @return The found solution for the transformation around this line
     */
    DualNumberAlgebra::DualNumber acos3(const UnitLine &a, const UnitLine &b) const noexcept;

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
