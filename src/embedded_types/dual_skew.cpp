//
// Created by sba on 19.07.21.
//

#include "dual_skew.h"
#include "unit_line.h"
#include "vector.h"

#include "precision.h"

DualSkew::DualSkew(const SkewMatrix &real, const SkewMatrix &dual) noexcept: DualEmbeddedMatrix(real, dual) {}

DualSkew::DualSkew(const UnitLine &line) noexcept: DualSkew(SkewMatrix(line.n()), SkewMatrix(line.m())) {}

UnitLine
DualSkew::screw() const noexcept {
    return UnitLine(
        UnitDirectionVector( // check if vector is a unit direction vector
            Vector( // create vector from skew matrix
                SkewMatrix( // check if matrix is skew matrix
                    Matrix3(this->data.topLeftCorner(3,3))
                )
            )
        ),
        MomentVector( // make vector a moment vector
            Vector( // create vector from skew matrix
                SkewMatrix( // check if matrix is skew matrix
                    Matrix3(this->data.bottomLeftCorner(3,3))
                )
            )
        )
    );
}

bool operator==(const DualSkew &lhs, const DualSkew &rhs) {
    return lhs.data.isApprox(rhs.data, Compare::instance().get_precision());
}
