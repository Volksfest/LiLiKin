//
// Created by sba on 19.07.21.
//

#ifndef DUAL_ALGEBRA_KINEMATICS_DUAL_EMBEDDED_MATRIX_H
#define DUAL_ALGEBRA_KINEMATICS_DUAL_EMBEDDED_MATRIX_H

#include <eigen3/Eigen/Eigen>

class Matrix3;
namespace DualNumberAlgebra {
    class DualNumber;
}

class RotationMatrix;
class Screw;

class DualFrame;
class DualSkewProduct;

/**
 * \brief A wrapper for an embedded 6x6 Matrix
 *
 * This type represents a dual number matrix as following:
 *  \f$
      \hat{E} = \begin{pmatrix}
       R & 0\\
       D & R
      \end{pmatrix}\f$
 * with \f$ E = R + \epsilon D \quad E \in \mathbb{D}^{3x3} \quad R,D \in \mathbb{R}^{3x3}\f$
 */
class DualEmbeddedMatrix {
public:
    /**
     * \brief Typedef for the Eigen Column Vector as 3x1 matrix
     */
    using Mat6 = Eigen::Matrix<double, 6, 6>;
protected:
    /**
     * \brief Internal data storage as Eigen type
     */
    Mat6 data;

    /**
     * \brief Protected constructor from Eigen type
     * Protected because the constructor does not have any type checks and an embedded matrix has some constraints.
     * See detailed type description
     * @param data Raw data in Eigen type
     */
    explicit DualEmbeddedMatrix(const Mat6 &data) noexcept;
public:
    /**
     * \brief Embed a dual representing 3x3 matrix to an other real representing 3x3 matrix within a 6x6 matrix
     *
     * See detailed type description
     *
     * @param real The real matrix
     * @param dual The dual matrix
     */
    DualEmbeddedMatrix(const Matrix3 &real, const Matrix3 &dual) noexcept;

    /**
     * \brief Represent a dual number as a 6x6 matrix
     *
     * This is done by by creating identity 3x3 matrices scaled by the real or dual part.
     * Then embed the identities
     *
     *  \f$
      \begin{pmatrix}
       r \cdot I_3 & 0\\
       d \cdot I_3 & r \cdot I_3
      \end{pmatrix}\f$
     * @param dn The represented dual number
     */
    explicit DualEmbeddedMatrix(const DualNumberAlgebra::DualNumber &dn) noexcept;

    /**
     * \brief The negation of the matrix
     * @return The negated matrix
     */
    DualEmbeddedMatrix operator-() const noexcept;

    /**
     * \brief The sum of the matrices
     * @param rhs The right-hand-side
     * @return The sum
     */
    DualEmbeddedMatrix operator+(const DualEmbeddedMatrix &rhs) const noexcept;

    /**
     * \brief The difference of the matrices
     * @param rhs The right-hand-side
     * @return The difference
     */
    DualEmbeddedMatrix operator-(const DualEmbeddedMatrix &rhs) const noexcept;

    /**
     * \brief The matrix multiplication
     * @param rhs The right-hand-side
     * @return The product
     */
    DualEmbeddedMatrix operator*(const DualEmbeddedMatrix &rhs) const noexcept;

    /**
     * \brief Return the corresponding Eigen data type (READ-ONLY)
     * @return The Eigen type of the matrix
     */
    const Mat6 & get() const noexcept;

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
};


#endif //DUAL_ALGEBRA_KINEMATICS_DUAL_EMBEDDED_MATRIX_H
