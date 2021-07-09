//
// Created by sba on 08.07.21.
//

#ifndef DAK_HOMOGENOUS_H
#define DAK_HOMOGENOUS_H

#include <eigen3/Eigen/Eigen>

class Vector;
class RotationMatrix;
class AdjungateMatrix;

//class HomogenousMatrix;
//std::ostream &operator<<(std::ostream &stream, const HomogenousMatrix &m);

class HomogenousMatrix {
private:
    using Mat4 = Eigen::Matrix<double,4,4>;

    Mat4 data;

    HomogenousMatrix(Mat4 data) noexcept;
public:
    HomogenousMatrix(const RotationMatrix & rot, const Vector & vector) noexcept;
    explicit HomogenousMatrix(const AdjungateMatrix &adj) noexcept;

    HomogenousMatrix operator*(const HomogenousMatrix &rhs) const noexcept;
    HomogenousMatrix inverse() const noexcept;

    RotationMatrix R() const noexcept;
    Vector p() const noexcept;

    friend std::ostream &operator<<(std::ostream &stream, const HomogenousMatrix &m);
};



#endif //DAK_HOMOGENOUS_H
