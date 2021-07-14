//
// Created by sba on 06.07.21.
//

#ifndef DAK_MATRIX6_H
#define DAK_MATRIX6_H

#include <eigen3/Eigen/Eigen>

#include "dual_number.h"

class Matrix3;
class Pluecker;
class RotationMatrix;
class HomogenousMatrix;

class Vector;
class AdjungateMatrix;
/**
* \brief A weaker wrapper for a generic 6x6 Matrix
*/
class Matrix6 {
private:
    using Mat6 = Eigen::Matrix<double, 6, 6>;

    Eigen::Matrix<double,6,6> data;

    explicit Matrix6(Mat6 data) noexcept;
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
    Matrix6(const Matrix3 &real, const Matrix3 &dual) noexcept;

    /**
     * \brief Represent a dual number as a 6x6 matrix
     *
     * This is done by by creating identity 3x3 matrices scaled by the real or dual part.
     * Then embed the dual part scaled identity into the real part scaled identitty.
     *
     *  \f$
      \begin{pmatrix}
       r \cdot I_3 & 0\\
       d \cdot I_3 & r \cdot I_3
      \end{pmatrix}\f$
     * @param dn The represented dual number
     */
    explicit Matrix6(const DualNumberAlgebra::DualNumber &dn) noexcept;

    /**
     * \brief The square operation as 6x6 matrix of a pluecker line
     * @param pl The pluecker line
     */
    explicit Matrix6(const Pluecker &pl) noexcept;

    /**
     * \brief The negation of the matrix
     * @return The negated matrix
     */
    Matrix6 operator-() const noexcept;

    /**
     * \brief The sum of the matrices
     * @param rhs The right-hand-side
     * @return The sum
     */
    Matrix6 operator+(const Matrix6 &rhs) const noexcept;

    /**
     * \brief The difference of the matrices
     * @param rhs The right-hand-side
     * @return The difference
     */
    Matrix6 operator-(const Matrix6 &rhs) const noexcept;

    /**
     * \brief The matrix multiplication
     * @param rhs The right-hand-side
     * @return The product
     */
    Matrix6 operator*(const Matrix6 &rhs) const noexcept;

//    /*
//     * \brief The matrix-vector-multiplication with a pluecker line as a map
//     *
//     * This actually only works with matrices representing dual numbers which is the case here.
//     *
//     * @param rhs To be mapped pluecker line
//     * @return The mapped pluecker line
//     */
//    Pluecker operator*(const Pluecker &rhs) const noexcept;

    friend Pluecker operator*(const Matrix6 &lhs, const Pluecker &rhs) noexcept;

    friend AdjungateMatrix;
};



/**
 * \brief Representation of a frame as a 6x6 matrix
 *
 * Can be interpreted as a dual number of a rotation matrix.
 * The dual part gives the translation but actually it is interlaced with the rotation
 */
class AdjungateMatrix : public Matrix6 {
private:
    explicit AdjungateMatrix(const Mat6 &mat) noexcept;
    explicit AdjungateMatrix(const Matrix6 &mat) noexcept;
public:
    /**
     * \brief Create a 6x6 matrix with the orientation and position
     *
     * \todo Change second parameter to Pointvector
     *
     * theoretically a dual number of rotation matrices
     * @param real The rotation matrix
     * @param dual The product of the skew translation matrix and the rotation matrix (pxR)
     */
    AdjungateMatrix(const RotationMatrix &real, const Matrix3 &dual) noexcept;

    /**
     * \brief Create the frame representation of a 6x6 matrix by a homogeneous 4x4 matrix
     * @param pl The homogeneous matrix to transform
     */
    explicit AdjungateMatrix(const HomogenousMatrix &pl) noexcept;

    /**
     * \brief Frame multiplication yields to another frame
     *
     * Frame A to ref W multiplied with Frame B to ref A yields to the Frame B to ref W
     * @param rhs The frame to transform
     * @return The transformed frame
     */
    AdjungateMatrix operator*(const AdjungateMatrix &rhs) const noexcept;

    /**
     * \brief The frame inversion
     * @return The inverted frame
     */
    AdjungateMatrix inverse() const noexcept;

    /**
     * \brief Retrieve the rotation matrix of the frame
     * @return The rotation matrix
     */
    RotationMatrix R() const noexcept;

    /**
     * \brief Retrieve the skewed translation rotation product
     * @return Return pxR
     */
    Matrix3 pxR() const noexcept;

    /**
     * \brief Retrieve the translation of the frame
     *
     * This is actually done by internally extracted it from pxR
     * @return The translation(position) of the frame
     */
    Vector p() const noexcept;

    friend Pluecker;
};

#endif //DAK_MATRIX6_H
