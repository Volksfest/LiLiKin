//
// Created by sba on 06.07.21.
//

#ifndef DAK_PLUECKER_H
#define DAK_PLUECKER_H

#include <eigen3/Eigen/Eigen>
#include "dual_number.h"

class Pluecker;
class SkewMatrix;
class Matrix3;
class Matrix6;
class AdjungateMatrix;

class PointVector;
class MomentVector;
class DirectionVector;
class HomogenousMatrix;

/**
 * \brief Wrapper class for Eigen 3x1 Matrix (which is a Vector)
 *
 * Just defined to have semantic data types.
 */

class Vector {
private:
    using Vec3 = Eigen::Matrix<double,3,1>;

    Vec3 data;

    explicit Vector(const Vec3 &data) noexcept;
public:
    /**
     * \brief Default constructor for zero vector
     */
    Vector() = default;

    /**
     * \brief Construct a vector with its three elements.
     *
     * @param a First element
     * @param b Second element
     * @param c Thrid element
     */
    Vector(double a, double b, double c) noexcept;

    /**
     * \brief Construct a Vector from a skewmatrix
     *
     * For a vector w a skewmatrix S[w] represets the lefthandside of the crossproduct:
     *
     * \f$
       w \: \times\: v = S[w] \cdot v
       \f$
     *
     *
       \f$
      S[w] =
      \begin{pmatrix}
       0 & -w_2 & w_1 \\
       w_2 & 0  & -w_0 \\
       -w_1 & w_0 & 0
       \end{pmatrix}\f$
     * @param skew Skewmatrix representation of the vector to create
     */
    Vector(const SkewMatrix &skew) noexcept;

    /**
     * \brief Operator overload for addition and assigning
     *
     * @param rhs Right-hand-side summand
     * @return The sum which is the object itself
     */
    Vector & operator+=(const Vector &rhs) noexcept;

    /**
     * \brief Operator overload for addition
     * @param rhs  Right-hand-side summand
     * @return The sum vector
     */
    Vector operator+(const Vector &rhs) const noexcept;

    /**
     * \brief Operator overload for substraction
     * @param rhs Right-hand-side operand (subtrahend)
     * @return The difference vector
     */
    Vector operator-(const Vector &rhs) const noexcept;

    /**
     * \brief Operator overload for scalar multiplication
     *
     * This is the linear scaling of a vector.
     * There is also a non-member function for the operands being exchanged to allow commutative operation.
     *
     * @param rhs Scaling factor
     * @return The scaled vector
     */
    Vector operator*(double rhs) const noexcept;

    /**
     * \brief Operator overload for scalar division
     *
     * This is also a linear scaling of a vector.
     *
     * \except std::logic_error If the right-hand-side is zero
     *
     * @param rhs Scaling factor as an inverse
     * @return The scaled vector
     */
    Vector operator/(double rhs) const;

    /**
     * \brief Operator overload for dot product
     *
     * This can be misinterpreted but is somehow also a usual notation.
     * @param rhs Second vector for dot product
     * @return The scaler dot product
     */
    double operator*(const Vector &rhs) const noexcept;

    /**
     * \brief The cross product of two vectors
     * @param rhs The second vector for the cross product
     * @return The cross product vector
     */
    Vector cross(const Vector &rhs) const noexcept;

    /**
     * \brief Calculate the norm of the vector
     * @return The norm of the vector
     */
    double norm() const noexcept;

    /**
     * \brief Cacluate the normalized vector
     *
     * \except std::logic_error If the norm is zero
     *
     * @return The normalized vector
     */
    Vector normal() const;

    /**
     * \brief Check if vector is the zero vector
     * @return True if zero vector
     */
    bool is_zero() const noexcept;

    friend Pluecker;
    friend PointVector;
    friend MomentVector;
    friend DirectionVector;
    friend Matrix3;
    friend SkewMatrix;
    friend HomogenousMatrix;
};

/**
 * \brief Helper function for commutative scaling of vector
 *
 * @see Vector::operator*
 *
 * @param lhs Left-hand-side as the factor
 * @param rhs right-hand-side as the vector
 * @return The scaled vector
 */
Vector operator*(double lhs, const Vector &rhs) noexcept;

/**
 * \brief Helper function for better readable cross operations
 * @param lhs Left-hand-side vector
 * @param rhs Right-hand-side vector
 * @return The cross product vector
 */
Vector cross(const Vector &lhs, const Vector &rhs) noexcept;

struct ProjectionTrio;

/**
 * \brief Semantic vector for defining a vector
 */
class PointVector : public Vector {
public:
    /**
     * \brief Conversion constructor from a Vector
     *
     * @param rhs Vector to treat as point vector
     */
    explicit PointVector(const Vector &rhs) noexcept;

    /**
     * \brief Create a point vector directly from its elements
     * @param a First element
     * @param b Second element
     * @param c Third element
     */
    PointVector(double a, double b, double c) noexcept;
};

/**
 * \brief Semantic vector for a direction
 *
 * The norm of this vector has to be one, thus direction vectors are always unit vectors
 */
class DirectionVector : public Vector {
public:
    /**
     * \brief Conversion constructor from a Vector
     * \except std::logic_error If the vector is not an unit vector
     * @param rhs Vector to treat as a direction vector
     */
    explicit DirectionVector(const Vector &rhs);

    /**
     * \brief Create a direction vector directly from its elements
     * \except std::logoic_error If the vector is not an unit vector
     * @param a First element
     * @param b Second element
     * @param c Third element
     */
    DirectionVector(double a, double b, double c);
};

/**
 * \brief Semantic vector for a moment
 */
class MomentVector : public Vector {
public:
    /**
     * \brief Conversion constructor from a Vector
     * @param rhs Vector to treat as a moment
     */
    explicit MomentVector(const Vector &rhs) noexcept;
    /**
     * \brief Create a moment vector directly from its elements
     * @param a First element
     * @param b Second element
     * @param c Third element
     */
    MomentVector(double a, double b, double c) noexcept;
};

/**
 * \brief The representation of a Pluecker line
 *
 * The line is define as a 6x1 Eigen Matrix.
 * The head three elements define the direction and the last three the moment.
 *
 * \todo Pluecker is badly defined right now. First of all a dual vector is better
 * \todo create a Unit Line subclass such that normalize makes sense
 * \todo reformulate as Generic Line class
 * \todo create Generic Twist Class such align makes sense (Twist as generic non-zero pitch twist) and make Generic line a subclass
 * \todo create Unit Twist subclass such that normalize and align makes sense
 */
class Pluecker{
private:
    using Vec6 = Eigen::Matrix<double,6,1>;

    Vec6 data;
    bool is_valid;

    explicit Pluecker(Vec6 data) noexcept;
    Pluecker() = default;
public:
    /**
     * \brief Create a line going throw the points a and b
     * @param a The point a
     * @param b The point b
     */
    explicit Pluecker(const PointVector &a, const PointVector &b) noexcept;

    /**
     * \brief Create a line with a direction and an anchor
     * @param n The direction vector (unit vector!)
     * @param a The anchor point
     */
    explicit Pluecker(const DirectionVector &n, const PointVector &a) noexcept;

    /**
     * \brief Create a line with a direction and a moment according to pluecker coordinates
     * \todo actually for a Line the moment vector needs to be orthogonal and has only 2 DoF. Like this it would yield to a non-zero pitch unit twist. Unit is given by the unit direction vector
     *
     * @param n The direction vector (unit vector!)
     * @param m The moment vector according to pluecker coordinates
     */
    explicit Pluecker(const DirectionVector &n, const MomentVector &m) noexcept;

    /**
     * \brief Return the direction vector
     * \todo change to DirectionVector
     * @return The direction vector
     */
    Vector n() const noexcept;

    /**
     * \brief Return the moment vector
     * \todo change to MomentVector
     * @return  The moment vector
     */
    Vector m() const noexcept;

    /**
     * \brief Calculate the scalar dual number product between lines as an operator overload
     *
     * \todo naming! change to a concrete name
     *
     * Scalar product is a dual number and indicates the distance and angle between the lines
     *
     * @param rhs
     * @return
     */
    DualNumberAlgebra::DualNumber operator*(const Pluecker &rhs) const noexcept;

    /**
     * \brief Retrieve the canoncical anchor of the line
     *
     * In short: retrieve an anchor.
     * Through the nature of the moment (for a unit and non-pitch line!) this is the point with the shortest distance to origin.
     * \todo create a parametrized function to get any point on a pluecker line
     * \todo create a designated anchor point for special class of line geometry with points!
     * \todo return a PointVector
     *
     * @return Canonical Anchor
     */
    Vector get_canonical_anchor() const noexcept;

    /**
     * \brief Get direction vector as a skew
     * \todo actually obsolet, let the user do it himself
     * @return skew matrix of the direction
     */
    SkewMatrix get_direction_skew() const noexcept;

    /**
     * \brief Calculate the transformatin (translation and rotation) around the line
     *
     * Using the generalized Rodriguez formula.
     *
     * @param angle The dual angle for translation and rotation amount
     * @return The adjungated frame which is the result of the transformation
     */
    AdjungateMatrix create_transform(DualNumberAlgebra::DualNumber angle) const noexcept;

    /**
     * \brief The negation and thus invert of direction of the line
     *
     * The lines are oriented as the direction vector can be negated.
     * Thus the lines are actually spears
     * \todo naming with lines and spears
     * @return The negated line
     */
    Pluecker operator-() const noexcept;

    /**
     * \brief this one makes actually no sense semantically but it is needed algebraically
     * \todo maybe private?
     * @param rhs
     * @return
     */
    Pluecker operator-(const Pluecker &rhs) const noexcept;

    /**
     * \brief Align the line thus remove the pitch of the screw
     *
     * This yields to a moment vector which is orthogonal to the lines direction
     *
     * @return Aligned line
     */
    Pluecker align() const;

    /**
     * \brief Normalize the line
     *
     * Probably this is unnecessary as the direction is already a unit and will always be!
     * @return
     */
    Pluecker normalize() const;

    /**
     * \brief Create a decomposition of a line as projection to this line
     * @param l Line to decompose
     * @return The decomposition in a decomposition
     */
    ProjectionTrio project(const Pluecker &l) const noexcept;

    /**
     * \brief Project the anchor of a line to the this line
     * \todo Probably change to PointVector
     * @param l The line where the anchor should be projected
     * @return The projection point
     */
    Vector point_project(const Pluecker &l) const noexcept;

    /**
     * \brief Create a translated line going through the new anchor
     * \todo Change to PointVector
     * @param new_anchor The new anchor
     * @return The new parallel line goint through the new anchor
     */
    Pluecker parallel_through_anchor(const Vector &new_anchor) const noexcept;

    friend Matrix6;
    friend AdjungateMatrix;
    friend Pluecker operator*(const Matrix6 &lhs, const Pluecker &rhs) noexcept;
};

/**
 * \brief A struct containing the projection of a line
 * \todo rename to LineDecomposition?
 * \todo add ref to reference line?
 * \todo change types. These lines are not unit non-zero pitch lines!
 */
struct ProjectionTrio {
Pluecker p; //!< projected line to the reference
Pluecker r; //!< rejected line orthogonal to reference but in the plane spanned by the line and the reference
Pluecker o; //!< orthogonal line orthogonal to reference and the line
};

#endif //DAK_PLUECKER_H
