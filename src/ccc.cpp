//
// Created by sba on 08.07.21.
//

#include "ccc.h"

using namespace DualNumberAlgebra;

CCCMechanism::CCCMechanism(const UnitLine &l12, const UnitLine &l23, const UnitLine &l34, const DualFrame &zero_posture) noexcept
    : l12(l12), l23(l23), l34(l34), zero_posture(zero_posture) {}

DualFrame
CCCMechanism::forward(const Configuration &config) const noexcept {
    return
        DualFrame(DualSkewProduct(this->l12, config.phi_1)) *
        DualFrame(DualSkewProduct(this->l23, config.phi_2)) *
        DualFrame(DualSkewProduct(this->l34, config.phi_3)) *
        this->zero_posture;
}

std::vector<Configuration>
CCCMechanism::inverse(const DualFrame &pose) const {

    auto s = pose * this->zero_posture.inverse();

    DualSkew crossterm(this->l23);
    DualEmbeddedMatrix uniterm = - crossterm * crossterm;
    DualEmbeddedMatrix squareterm = DualEmbeddedMatrix(1) - uniterm;

    DualNumber a = this->l12 * (uniterm * this->l34);
    DualNumber b = this->l12 * (crossterm * this->l34);
    DualNumber c = this->l12 * ((s - squareterm) * this->l34);

    auto phi2_solutions = solve_trigonometric_equation(a, b, c);
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
        DualFrame m2(DualSkewProduct(this->l23, phi_2));

        auto phi_1 = (-this->l12).acos3(s * this->l34, m2 * this->l34);
        auto phi_3 = this->l34.acos3(s.inverse() * this->l12, m2.inverse() * this->l12);

        // TODO change hack
        if (enable_offset) {
            auto intermediate_line = DualFrame(DualSkewProduct(this->l12,phi_1)) * m2 * this->l23;
            auto offset = l34.acos3(
                    intermediate_line,
                    intermediate_line.parallel_through_anchor(s.p()));
            phi_3 = phi_3 + offset;
        }

        solutions.push_back({phi_1, phi_2, phi_3});
    }

    return solutions;
}