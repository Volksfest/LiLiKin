//
// Created by sba on 19.07.21.
//

#ifndef DUAL_ALGEBRA_KINEMATICS_DUAL_FRAME_H
#define DUAL_ALGEBRA_KINEMATICS_DUAL_FRAME_H

#include "embedded_types/dual_embedded_matrix.h"
#include "embedded_types/dual_skew_product.h"
#include "base/vector.h"

class DualFrame: public DualEmbeddedMatrix {
private:
    explicit DualFrame(const DualEmbeddedMatrix &mat) noexcept;
protected:
    explicit DualFrame(const Mat6 &mat) noexcept;
public:
    DualFrame(const RotationMatrix &real, const PointVector &dual) noexcept;

    DualFrame(const DualSkewProduct &skew) noexcept;

    DualFrame operator*(const DualFrame &rhs) const noexcept;

    /**
     * \brief The frame inversion
     * @return The inverted frame
     */
    DualFrame inverse() const noexcept;

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

    friend UnitScrew operator*(const DualFrame &lhs, const UnitScrew &rhs) noexcept;
    friend Line operator*(const DualFrame &lhs, const Line &rhs) noexcept;
    friend UnitLine operator*(const DualFrame &lhs, const UnitLine &rhs) noexcept;

    friend std::ostream &operator<<(std::ostream &stream, const DualFrame &d);
};




#endif //DUAL_ALGEBRA_KINEMATICS_DUAL_FRAME_H
