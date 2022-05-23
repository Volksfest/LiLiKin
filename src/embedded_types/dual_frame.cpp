//
// Created by sba on 19.07.21.
//

#include "dual_frame.h"

#include "precision.h"

#include <iomanip>

DualFrame::DualFrame(const Mat6 &mat) noexcept: DualEmbeddedMatrix(mat) {}

DualFrame::DualFrame(const RotationMatrix &rot, const PointVector &trans) noexcept:
    DualEmbeddedMatrix(rot, SkewMatrix(trans) * rot) {}

DualFrame::DualFrame(const DualSkewProduct &argument) noexcept: DualEmbeddedMatrix(argument.skew()) {
    auto orthoterm = argument.skew();
    auto uniterm = - orthoterm * orthoterm;
    auto nullterm = DualEmbeddedMatrix(1.0) - uniterm;
    auto angle = argument.angle();

    this->data = (DualEmbeddedMatrix(cos(angle)) * uniterm +
                           DualEmbeddedMatrix(sin(angle)) * orthoterm +
                           nullterm).get();
}

DualFrame DualFrame::operator*(const DualFrame &rhs) const noexcept {
    return DualFrame(this->data * rhs.data);
}

DualFrame DualFrame::inverse() const noexcept {
    Eigen::Matrix<double, 3,3> rot = this->R().get().transpose();
    Eigen::Matrix<double, 3,3> dual = this->pxR().get().transpose();

    Mat6 data;
    data.topLeftCorner(3, 3)     = rot;
    data.topRightCorner(3, 3)     = Eigen::Matrix<double,3,3>::Zero(3, 3);
    data.bottomLeftCorner(3, 3)     = dual;
    data.bottomRightCorner(3, 3)     = rot;

    return DualFrame(data);
}

RotationMatrix DualFrame::R() const noexcept {
    return RotationMatrix(this->data.topLeftCorner(3,3));
}

Matrix3 DualFrame::pxR() const noexcept {
    return Matrix3(this->data.bottomLeftCorner(3,3));
}

PointVector DualFrame::p() const noexcept {
    return PointVector(Vector(SkewMatrix(this->pxR() * this->R().inverse())));
}

std::ostream &operator<<(std::ostream &stream, const DualFrame &d) {
    auto R = d.data.topLeftCorner(3,3);
    auto ps = d.data.bottomLeftCorner(3,3) * R.inverse();
    Eigen::Matrix<double,3,1> p;
    p << ps(2,1), ps(0,2), ps(1,0);
    stream << std::fixed;
    stream << "┌─────────────────────────────────┬───────────┐" << std::endl;
    for( int i = 0; i < 3; i++) {
        stream << "│";
        for (int j = 0; j < 3; j++) {
            stream << std::setw(10) << R(i,j) << " ";
        }
        stream << "│" << std::setw(10) << p(i) << " │" << std::endl;
    }
    stream << "└─────────────────────────────────┴───────────┘" << std::endl;
    stream << std::defaultfloat;
    return stream;
}

bool operator==(const DualFrame &lhs, const DualFrame &rhs) noexcept {
    auto check_R = (lhs.R().get() - rhs.R().get()).array().abs();
    if ( ! (check_R< Compare::instance().get_precision()).all() ) {
        return false;
    }
    auto check_pxR = (lhs.pxR().get() - rhs.pxR().get()).array().abs();
    if ( ! (check_pxR< Compare::instance().get_precision()).all() ) {
        return false;
    }
    return true;
}

DualSkewProduct DualFrame::constructive_line() const noexcept {
    auto diff = *this - this->inverse();
    auto n_skew = diff.get().topLeftCorner(3,3);
    auto m_skew = diff.get().bottomLeftCorner(3,3);

    auto R = this->data.topLeftCorner(3,3);
    auto pxR = this->data.bottomLeftCorner(3,3);

    auto screw = Screw(
            DirectionVector(Vector(SkewMatrix(Matrix3(n_skew)))),
            MomentVector(Vector(SkewMatrix(Matrix3(m_skew)))));

    auto angle = acos(0.5 * (DualNumberAlgebra::DualNumber(R.trace(), pxR.trace()) - 1));

    return {screw.to_line(), angle};
}