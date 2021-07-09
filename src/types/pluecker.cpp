//
// Created by sba on 06.07.21.
//

#include "types.h"


Pluecker::Pluecker(Vec6 data) noexcept {
    this->data = data;
    this->is_valid = (data.head(3).norm() >0.5);
}

Pluecker::Pluecker(const PointVector &a, const PointVector &b) noexcept {
    Vector n = b - a;
    if (n.is_zero()) {
        this->is_valid = false;
    } else {
        n = n / n.norm();
        Vector m = cross(a, n);
        this->data << n.data, m.data;
        this->is_valid = true;
    }
}

Pluecker::Pluecker(const DirectionVector &n, const PointVector &a) noexcept {
    if (n.is_zero()) {
        this->is_valid = false;
    } else {
        Vector m = cross(a, n);
        this->data << n.data, m.data;
        this->is_valid = true;
    }
}

Pluecker::Pluecker(const DirectionVector &n, const MomentVector &m) noexcept {
    if (n.is_zero()) {
        this->is_valid = false;
    } else {
        this->data << n.data, m.data;
        this->is_valid = true;
    }
}

Vector
Pluecker::n() const noexcept {
    return Vector(this->data.head(3));
}

Vector
Pluecker::m() const noexcept {
    return Vector(this->data.tail(3));
}

DualNumberAlgebra::DualNumber
Pluecker::operator*(const Pluecker &rhs) const noexcept {
    using DualNumberAlgebra::DualNumber;
    Eigen::Matrix<DualNumber,3,1> l, r;

    l<< DualNumber(this->data(0,0), this->data(3,0)),
            DualNumber(this->data(1,0), this->data(4,0)),
            DualNumber(this->data(2,0), this->data(5,0));

    r<< DualNumber(rhs.data(0,0), rhs.data(3,0)),
            DualNumber(rhs.data(1,0), rhs.data(4,0)),
            DualNumber(rhs.data(2,0), rhs.data(5,0));

    return l.transpose() * r;
}

Vector
Pluecker::get_canonical_anchor() const noexcept {
    return cross(this->n(), this->m());
}

SkewMatrix
Pluecker::get_direction_skew() const noexcept {
    return SkewMatrix(this->n());
}

AdjungateMatrix
Pluecker::create_transform(DualNumberAlgebra::DualNumber angle) const noexcept {
    Matrix6 orthoterm(*this);
    Matrix6 uniterm = - orthoterm * orthoterm; //eq -i^2
    Matrix6 nullterm = Matrix6(1.0) - uniterm;

    return AdjungateMatrix(Matrix6(cos(angle)) * uniterm +
                           Matrix6(sin(angle)) * orthoterm +
                           nullterm);
}

Pluecker
Pluecker::align() const {
    if (!this->is_valid) {
        throw std::logic_error("This is not a valid line");
    }
    Pluecker a;

    Eigen::Matrix<double,3,1> n,m;
    n = this->data.head(3);
    m = this->data.tail(3);

    a.data.head(3) = n;
    a.data.tail(3) = m - (n.dot(m) / n.dot(n)) * n;
    return a;
}

Pluecker
Pluecker::normalize() const {
    if (!this->is_valid) {
        throw std::logic_error("This is not a valid line");
    }

    double norm = this->n().norm();

    return Pluecker(this->data / norm);
}

ProjectionTrio
Pluecker::project(const Pluecker &l) const noexcept {

    Matrix6 adj(*this);

    Pluecker p, r, o;

    o = adj * l;
    r = -adj * o;
    p = l - r;

    if (!o.is_valid) {
        auto anchor_n = this->get_canonical_anchor();
        auto c = anchor_n - l.get_canonical_anchor();

        auto n = cross(cross(this->n(), c), this->n());

        o = Pluecker(DirectionVector(n.normal()),
                     PointVector(anchor_n/n.norm()));
    }

    return {p, r, o};
}

Vector
Pluecker::point_project(const Pluecker &l) const noexcept {
    Vector point_projection = this->get_canonical_anchor();
    double cross_norm = cross(l.n(), this->n()).norm();
    if (cross_norm != 0) {
        auto projection = l.get_direction_skew();
        point_projection += this->n() * (projection * projection * this->n() * (this->get_canonical_anchor() - l.get_canonical_anchor()) /
                                      cross_norm / cross_norm);
    }
    return point_projection;
}

// not really line geometry
Pluecker
Pluecker::parallel_through_anchor(const Vector &new_anchor) const noexcept {
    return Pluecker(
            DirectionVector(this->n()),
            PointVector(new_anchor)
            );
}

Pluecker
Pluecker::operator-(const Pluecker &rhs) const noexcept {
    return Pluecker(this->data - rhs.data);
}

Pluecker
Pluecker::operator-() const noexcept {
    return Pluecker(-this->data);
}
