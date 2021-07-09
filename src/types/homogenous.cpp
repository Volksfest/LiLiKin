//
// Created by sba on 08.07.21.
//

#include "types/homogenous.h"
#include "types/pluecker.h"
#include "types/matrix3.h"
#include "types/matrix6.h"

HomogenousMatrix::HomogenousMatrix(Mat4 data) noexcept {
    this->data = data;
}


HomogenousMatrix::HomogenousMatrix(const RotationMatrix & rot, const Vector & vector) noexcept {
    this->data = Mat4::Identity(4,4);
    this->data.topLeftCorner(3,3) = rot.data;
    this->data.topRightCorner(3,1) = vector.data;
}

HomogenousMatrix::HomogenousMatrix(const AdjungateMatrix &adj) noexcept
    : HomogenousMatrix(adj.R(), adj.p() ) {}

HomogenousMatrix
HomogenousMatrix::operator*(const HomogenousMatrix &rhs) const noexcept {
    return HomogenousMatrix(this->data * rhs.data);
}

HomogenousMatrix
HomogenousMatrix::inverse() const noexcept {
    return HomogenousMatrix(this->data.inverse()); //TODO could be done faster
}

RotationMatrix
HomogenousMatrix::R() const noexcept {
    return RotationMatrix(this->data.topLeftCorner(3,3));
}

Vector
HomogenousMatrix::p() const noexcept {
    return Vector(this->data.topRightCorner(3,1));
}

std::ostream &operator<<(std::ostream &stream, const HomogenousMatrix &m) {
    return stream << m.data << std::endl;
}