//
// Created by sba on 08.07.21.
//

#ifndef DAK_HOMOGENOUS_H
#define DAK_HOMOGENOUS_H

#include <eigen3/Eigen/Eigen>


class Vector;
class RotationMatrix;
class AdjungateMatrix;

/**
 * \brief Homogenous 4x4 Matrix
 *
 * The implementation is a plain Eigen3 4x4 matrix.
 * No additional optimization is implemented
 */
class HomogenousMatrix {
private:
    using Mat4 = Eigen::Matrix<double,4,4>;

    Mat4 data;

    HomogenousMatrix(Mat4 data) noexcept;
public:
    /**
     * This creates a frame as homogenous matrix as:
     *
     * \f$
      F =
      \begin{pmatrix}
       R & t\\
       0 & a
      \end{pmatrix}\f$
     * @param[in] rot The rotation matrix of the frame
     * @param[in] vector The translation of the frame
     */
    HomogenousMatrix(const RotationMatrix & rot, const Vector & vector) noexcept;

    /**
     * This creates a frame as homogeneous matrix by converting it from the adjungated 6x6 matrix representation.
     *
     * Internally R and t are retrieved from the adjungated matrix and used to create the 4x4 matrix.
     * @param[in] adj The adjungate of the frame
     */
    explicit HomogenousMatrix(const AdjungateMatrix &adj) noexcept;

    /**
     * The multiplication between two homogeneous matrices
     *
     * The homogeneous matrices are a non commutative group within the multiplication.
     * Thus the result is also an homogeneous matrix.
     *
     * @param[in] rhs The right-hand-side of the multiplication.
     * @return The compound frame
     */
    HomogenousMatrix operator*(const HomogenousMatrix &rhs) const noexcept;

    /**
     * The inverse of the homogeneous matrix.
     *
     * This is not optimized and just the inverse in the Eigen matrix.
     *
     * @return The inverse of the homogeneous matrix.
     */
    HomogenousMatrix inverse() const noexcept;

    /**
     * The rotation part of the homogeneous matrix.
     *
     * @return The rotation part of the homogeneous matrix.
     */
    RotationMatrix R() const noexcept;

    /**
     * The translation part of the homogeneous matrix
     *
     * @return The translation part of the homogeneous matrix
     */
    Vector p() const noexcept;

    /**
     * The friend stream operation to display the homogeneous matrix in the standard output
     *
     * @param stream The stream to write to
     * @param[in] m The homogeneous matrix to be written
     * @return The stream again to chain the operation
     */
    friend std::ostream &operator<<(std::ostream &stream, const HomogenousMatrix &m);
};

#endif //DAK_HOMOGENOUS_H
