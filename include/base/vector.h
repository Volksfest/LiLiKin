//
// Created by sba on 19.07.21.
//

#ifndef DUAL_ALGEBRA_KINEMATICS_VECTOR_H
#define DUAL_ALGEBRA_KINEMATICS_VECTOR_H

#include <eigen3/Eigen/Eigen>

class Matrix3;
class SkewMatrix;
class Screw;

/**
 * \brief Wrapper class for Eigen 3x1 Matrix (which is a Vector)
 *
 * Just defined to have semantic data types.
 */
class Vector {
protected:
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
    explicit Vector(const SkewMatrix &skew) noexcept;

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

    friend Matrix3;
    friend SkewMatrix;
    friend Screw;
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

class UnitDirectionVector;

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

    UnitDirectionVector normal() const;
};

class UnitDirectionVector: public DirectionVector {
public:
    explicit UnitDirectionVector(const Vector &rhs);

    UnitDirectionVector(double a, double b, double c);

    UnitDirectionVector normal() const;
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

#endif //DUAL_ALGEBRA_KINEMATICS_VECTOR_H
