//
// Created by sba on 19.07.21.
//

#ifndef DUAL_ALGEBRA_KINEMATICS_DUAL_FRAME_H
#define DUAL_ALGEBRA_KINEMATICS_DUAL_FRAME_H

#include "embedded_types/dual_embedded_matrix.h"
#include "embedded_types/dual_skew_product.h"
#include "base/vector.h"

/**
 * \brief Frame representation as a dual embedded matrix
 *
 * The rotation matrix is the real part.
 * The skew translation with rotation is the dual part
 */
class DualFrame: public DualEmbeddedMatrix {
protected:
    /**
     * \brief Protected constructor with raw Eigen type
     *
     * It is protected to enforce a high level user to use the other constructor.
     * Especially important to avoid wrong rotation matrix format.
     *
     * @param mat Data in Eigen type
     */
    explicit DualFrame(const Mat6 &mat) noexcept;
public:
    /**
     * \brief Embed a rotation and a position to a dual frame
     *
     * @param real The real matrix representing the rotation
     * @param dual The dual matrix representing a mix of rotation and translation
     */
    DualFrame(const RotationMatrix &real, const PointVector &dual) noexcept;

    /**
     * \brief Map a DualSkewProduct to a frame
     *
     * This performs a generalized Rodriguez-Formula.
     * It corresponds to a map from a lie-algebra to a lie-group, allthough we are talking about se(3) and SE(3)
     *
     * @param skew The skewproduct to map
     */
    DualFrame(const DualSkewProduct &skew) noexcept;

    /**
     * \brief Multiplication/Concatenation of frames
     *
     * Dual frames form a group just like usual frames
     *
     * @param rhs The other frame
     * @return The concatenated frame
     */
    DualFrame operator*(const DualFrame &rhs) const noexcept;

    /**
     * \brief The frame inversion
     * @return The inverted frame
     */
    DualFrame inverse() const noexcept;

    /**
     * \brief Create the constructive line with its angle for the frame
     *
     * It is the inversion of the exponential map, thus the logarithm to calculate the lie algebra from the lie group.
     * So, the se(3) from the SE(3) is calculated.
     *
     * @return The transforming line with the angle yielding to this frame
     */
    DualSkewProduct constructive_line() const noexcept;

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
    PointVector p() const noexcept;

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

    /**
     * \brief Output formatting to homogenous matrices for a better readability
     * @param stream Output stream
     * @param d The frame to print
     * @return Output stream
     */
    friend std::ostream &operator<<(std::ostream &stream, const DualFrame &d);

    /**
     * \brief Comparison between two frames
     * @param lhs left-hand-side
     * @param rhs right-hand-side
     * @return True if equal
     */
    friend bool operator==(const DualFrame &lhs, const DualFrame &rhs) noexcept;
};

#endif //DUAL_ALGEBRA_KINEMATICS_DUAL_FRAME_H
