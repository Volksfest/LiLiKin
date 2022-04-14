//
// Created by sba on 06.07.21.
//

#ifndef DAK_MATRIX3_H
#define DAK_MATRIX3_H

#include <eigen3/Eigen/Eigen>

class Vector;
class DualEmbeddedMatrix;
class DualFrame;
class DualSkewProduct;

/**
 * \brief A wrapper for a generic 3x3 Matrix
 *
 */
class Matrix3 {
public:
    /**
     * Typedef for the Eigen 3x3 matrix
     */
    using Mat3 = Eigen::Matrix<double,3,3>;
protected:
    /**
     * Internal data as Eigen type
     */
    Eigen::Matrix<double,3,3> data;
public:
    /**
     * \brief Create the Matrix by a Eigen matrix
     * @param data Eigen type
     */
    explicit Matrix3(const Mat3 &data) noexcept;

    /**
     * Generic scalar matrix multiplication
     * @param rhs right factor
     * @return Scaled matrix
     */
    Matrix3 operator*(double rhs) const noexcept;

    /**
     * Generic matrix multiplication
     * @param rhs right factor
     * @return Product matrix
     */
    Matrix3 operator*(const Matrix3 &rhs) const noexcept;

    /**
     * Generic matrix-vector-multiplication
     *
     * Seeing the matrix as a linear map this maps the vector to another vector.
     * @param rhs Input vector
     * @return Mapped vector
     */
    Vector operator*(const Vector &rhs) const noexcept;

    /**
     * \brief Return the corresponding Eigen data type (READ-ONLY)
     * @return The Eigen type of the matrix
     */
    const Mat3 & get() const noexcept;
};

/**
 * Leftwise scalar multiplication
 *
 * Matrix scaling is commutative.
 * Thus a left-hand-side scalar version is also delivered.
 * @param lhs Scaling factor
 * @param rhs Input matrix
 * @return Scaled matrix
 */
Matrix3 operator*(double lhs, const Matrix3 &rhs) noexcept;

/**
 * \brief Orthonormal matrix representing rotations
 *
 * Right now only z-y-x convention is implemented to generate the matrix.
 *
 * Other possibility is to create a rotationmatrix via a Eigen3.
 * The matrix will be checked for being an orthogonal matrix.
 */
class RotationMatrix : public Matrix3 {
private:
    /**
     * \brief Protected constructor with raw Eigen type
     *
     * It is protected to enforce a high level user to use the other constructor to guarantee orthogonality
     *
     * @param data Data in Eigen type
     */
    explicit RotationMatrix(const Mat3 &data) noexcept;
public:
    /**
     * \brief Create a rotation matrix with the z-y-x convention.
     *
     * This is equivalent to Yaw(z) - Pitch(y) - Roll(x)
     * @param z Yaw angle
     * @param y Pitch angle
     * @param x Roll angle
     */
    RotationMatrix(double z, double y, double x) noexcept;

    /**
     * \brief Create a rotation matrix of raw data
     *
     * Matrix will be checked if it is a orthogonal matrix
     *
     * \exception std::domain_error Will be thrown if the input matrix is not an orthogonal matrix
     *
     * @param rhs The rotationmatrix
     *
     */
    static RotationMatrix RotationFromEigen(const Mat3 &data);

    /**
     * \brief Matrix-Matrix multiplication within rotation matrices.
     *
     * The set of ortho normal matrices are a non commutative group within multiplication.
     * Thus the product is also a rotation matrix
     * @param rhs The other rotation
     * @return The matrix representing the composed rotation
     */
    RotationMatrix operator*(const RotationMatrix &rhs) const noexcept;

    /**
     * \brief The inverse rotation.
     *
     * As the multiplication is a group, the inverse exists and is also within the rotation matrices.
     * Actually, it is easy to calculate as it is just the transpose
     * @return The inverse rotation
     */
    RotationMatrix inverse() const noexcept;

    /**
     * \brief A friend so that the DualFrame can create a RotationMatrix with Eigen types
     */
    friend DualFrame;
};

/**
 * \brief 3x3 matrices representing the skew of a vector
 *
 * This is primarily an Eigen3 3x3 matrix which gets a typecheck in its creation to ensure the skew constraint.
 *
 * The skew matrices are a commutative (abelian) group in addition but it is not implemented as it is not needed.
 * This is also missing in the super class.
 *
 * The same is also valid for the scaling.
 * But this is implemented in the parent thus it is possible to scale the skew matrix although the result sadly wouldn't be a SkewMatrix anymore.
 */
class SkewMatrix : public Matrix3 {
public:
    /**
     * Check if the input is a skew matrix and use it as skew matrix
     *
     * \exception std::domain_error Will be thrown if the input matrix is not a skew matrix
     *
     * @param rhs The matrix to be checked
     *
     */
    explicit SkewMatrix(const Matrix3 & rhs);

    /**
     * Create a skew matrix from a vector.
     *
     * A vector mapping - using the parents Matrix3::operator* - represents a cross product.
     * @param vector The input vector
     */
    explicit SkewMatrix(const Vector & vector) noexcept;
};

#endif //DAK_MATRIX3_H
