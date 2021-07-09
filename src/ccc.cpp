//
// Created by sba on 08.07.21.
//

#include "ccc.h"
#include "types/matrix6.h"
#include "dual_number.h"

#include "trigonometry_helper.h"

#include <iostream>

using namespace DualNumberAlgebra;

CCCMechanism::CCCMechanism(const Pluecker &l12, const Pluecker &l23, const Pluecker &l34, const AdjungateMatrix &zero_posture) noexcept
    : l12(l12), l23(l23), l34(l34), zero_posture(zero_posture) {}

AdjungateMatrix
CCCMechanism::forward(const Configuration &config) const noexcept {
    return
        this->l12.create_transform(config.phi_1) *
        this->l23.create_transform(config.phi_2) *
        this->l34.create_transform(config.phi_3) *
        this->zero_posture;
}

std::vector<Configuration>
CCCMechanism::inverse(const AdjungateMatrix &pose) const {

    AdjungateMatrix s = pose * this->zero_posture.inverse();

    Matrix6 crossterm(this->l23);
    Matrix6 uniterm = - crossterm * crossterm;
    Matrix6 squareterm = Matrix6(1) - uniterm;

    DualNumber a = this->l12 * (uniterm * this->l34);
    DualNumber b = this->l12 * (crossterm * this->l34);
    DualNumber c = this->l12 * ((s - squareterm) * this->l34);

    auto phi2_solutions = solve_trig(a, b, c);
    bool enable_offset = false;

    if (phi2_solutions.empty()) {
        auto p = s * this->l34;
        auto d = cross(this->l12.n(), this->l12.m() - p.m()).norm();

        phi2_solutions.emplace_back(DualNumber(0,-d));
        phi2_solutions.emplace_back(DualNumber(0,d));

        enable_offset = true;
    }

    std::vector<Configuration> solutions;

    for (auto phi_2 : phi2_solutions) {
        auto m2 = this->l23.create_transform(phi_2);

        auto phi_1 = acos3(s * this->l34, m2 * this->l34, -this->l12);
        auto phi_3 = acos3(s.inverse() * this->l12, m2.inverse() * this->l12, this->l34);

        // TODO change hack
        if (enable_offset) {
            auto intermediate_line = this->l12.create_transform(phi_1) * m2 * this->l23;
            auto offset = acos3(
                    intermediate_line,
                    intermediate_line.parallel_through_anchor(s.p()),
                    l34);
            phi_3 = phi_3 + offset;
        }

        solutions.push_back({phi_1, phi_2, phi_3});
    }

    return solutions;
}