//
// Created by sba on 22.09.20.
//

#ifndef DAK_MATH_H
#define DAK_MATH_H

using DualNumber = DualNumberAlgebra::DualNumber<double>;
using Pluecker = RobotAlgebra::Pluecker;
using Vector = RobotAlgebra::Vector;
using AdjungateMatrix = RobotAlgebra::AdjungateMatrix;
using HomMatrix = RobotAlgebra::HomogenousMatrix;

using DualVector = Eigen::Matrix<DualNumber,3,1>;

DualVector toDualLine(const Pluecker &p) noexcept;

DualNumber dot(const Pluecker &lhs, const Pluecker &rhs) noexcept;

DualNumber acos3(const Pluecker &a, const Pluecker &b, const Pluecker &n) noexcept;

bool solve(
        const AdjungateMatrix &ik,
        const AdjungateMatrix &zp,
        const Pluecker &l12,
        const Pluecker &l23,
        const Pluecker &l34,
        DualNumber &phi_1,
        DualNumber &phi_2,
        DualNumber &phi_3, bool other_solution) noexcept;

#endif //DAK_MATH_H

