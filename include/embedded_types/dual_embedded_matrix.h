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
* \brief A weaker wrapper for a generic 6x6 Matrix
*/
class DualEmbeddedMatrix {
public:
    using Mat6 = Eigen::Matrix<double, 6, 6>;
protected:
    Mat6 data;

    explicit DualEmbeddedMatrix(const Mat6 &data) noexcept;
public:
    /**
     * \brief Embed a dual representing 3x3 matrix to an other real representing 3x3 matrix within a 6x6 matrix
     *
     * This is done by using the upper left and lower right for the real part and the bottom left for the embedded dual part.
     * The upper right is zero.
     *
     *  \f$
      \begin{pmatrix}
       R & 0\\
       D & R
      \end{pmatrix}\f$
     *
     * @param real The real matrix
     * @param dual The dual matrix
     */
    DualEmbeddedMatrix(const Matrix3 &real, const Matrix3 &dual) noexcept;

    /**
     * \brief Represent a dual number as a 6x6 matrix
     *
     * This is done by by creating identity 3x3 matrices scaled by the real or dual part.
     * Then embed the dual part scaled identity into the real part scaled identity.
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

    const Mat6 & get() const noexcept;

    friend Screw operator*(const DualEmbeddedMatrix &lhs, const Screw &rhs) noexcept;
};


#endif //DUAL_ALGEBRA_KINEMATICS_DUAL_EMBEDDED_MATRIX_H
