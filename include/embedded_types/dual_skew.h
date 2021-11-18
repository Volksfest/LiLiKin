//
// Created by sba on 19.07.21.
//

#ifndef DUAL_ALGEBRA_KINEMATICS_DUAL_SKEW_H
#define DUAL_ALGEBRA_KINEMATICS_DUAL_SKEW_H

#include "embedded_types/dual_embedded_matrix.h"
#include "screws/unit_line.h"
#include "base/matrix3.h"

class DualSkewProduct;

/**
 * \brief The dual skew of an UnitLine
 *
 * It behaves similar to the SkewMatrix but for lines
 */
class DualSkew: public DualEmbeddedMatrix {
private:
    /**
     * \brief Generate the DualSkew by the real and dual skews
     *
     * Mostly given by a skew direction and a skew moment
     *
     * @param real The real part of the dual skew
     * @param dual The dual part of the dual skew
     */
    DualSkew(const SkewMatrix &real, const SkewMatrix &dual) noexcept;
public:
    /**
     * \brief create a DualSkew for a line
     *
     * Meaning a skew direction as real part and a skew moment for the dual part
     * @param line The line to transform to a skew
     */
    explicit DualSkew(const UnitLine &line) noexcept;

    /**
     * \brief Invert the skew to line again
     * @return The corresponding line
     */
    UnitLine screw() const noexcept;

    /**
     * \brief Comparison between two skews
     * @param lhs left-hand-side
     * @param rhs right-hand-side
     * @return True if equal
     */
    friend bool operator==(const DualSkew &lhs, const DualSkew &rhs);
};

#endif //DUAL_ALGEBRA_KINEMATICS_DUAL_SKEW_H
