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

std::tuple<DualFrame, UnitLine, UnitLine>
CCCMechanism::forward_verbose(const Configuration &config) const {
    auto fk1 = DualFrame(DualSkewProduct(this->l12, config.phi_1));
    auto fk12 = fk1 * DualFrame(DualSkewProduct(this->l23, config.phi_2));
    auto fk123 = fk12 * DualFrame(DualSkewProduct(this->l34, config.phi_3));

    return std::make_tuple(fk123 * this->zero_posture, fk1 * this->l23, fk12 * l34);
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

    DualNumber distance = this->l12.get_distance(s * this->l34);
    std::vector<DualNumber> phi2_solutions = solve_trigonometric_equation(a, b, c);

    std::vector<Configuration> solutions;

    for (auto phi_2 : phi2_solutions) {
        // coincide
        if ( (abs(distance.real()) < 0.000001 && distance.dual() < 0.000001) || (abs(distance.real()-M_PI) < 0.00001 && abs(distance.dual()) < 0.00001)) {

            if (abs(distance.real()-M_PI) < 0.0001) {
                phi_2 = phi_2 + M_PI;
            }
            DirectionVector u(cross(this->l12.n(), this->l23.n()));
            UnitLine orthogonal(
                u.normal(),
                PointVector(0,0,0)
            );

            DualFrame m2(DualSkewProduct(this->l23, phi_2));
            auto phi_1 = (-this->l12).acos3(s * orthogonal, m2 * orthogonal);
            auto phi_3 = DualNumber();

            solutions.push_back({phi_1, phi_2, phi_3});
        } else {
            // parallel
            if ( abs(distance.real()) < 0.000001 || abs(distance.real()-M_PI) < 0.000001) {
                DualNumber flip;
                if ( abs(distance.real()-M_PI) < 0.00001) {
                    flip = M_PI;
                }
                DualFrame m2(DualSkewProduct(this->l23, phi_2));
                DualFrame m2_rot(DualSkewProduct(this->l23, phi_2.real()));
                auto phi_1 = (-this->l12).acos3(s * this->l34, m2 * this->l34);
                auto phi_3 = this->l34.acos3(s.inverse() * this->l12, m2.inverse() * this->l12);

                DualFrame m1_combined(DualSkewProduct(this->l12,(phi_1 + phi_3).real()));
                auto pre_pose = m1_combined * m2_rot * this->zero_posture;

                // Some kind of Paden-Kahan Problem
                PointVector anchor = this->l12.get_canonical_anchor();
                PointVector p_target = this->l12.orthogonal_plane_projection(anchor,pose.p());
                PointVector p_begin = this->l12.orthogonal_plane_projection(anchor, (pre_pose).p());

                PointVector o_o = this->l12.get_canonical_anchor();
                PointVector o_r = (m1_combined * m2_rot * this->l34).point_project(o_o);

                Vector p = p_target - p_begin;
                Vector o = o_r - o_o;
                Vector n = (m1_combined * this->l23).n();
                Vector w = this->l12.n();

                double minus_p_2 = -(n * o);
                double q = - (p * ( p + 2 * o));
                double sq = sqrt(minus_p_2 * minus_p_2 - q);
                std::vector<double> alphas = {minus_p_2 - sq, minus_p_2 + sq};

                Vector v = p + o;
                for (double alpha : alphas) {
                    Vector u = o + alpha * n;
                    double phi_off = atan2(w * cross(u,v), u * v);
                    auto tmp_translation = this->l12.point_project(pose.p()) - this->l12.point_project((m2 * this->zero_posture).p());
                    double translation_off = tmp_translation.norm() * (tmp_translation * this->l12.n() > 0 ? 1.0 : -1.0);

                    DualFrame m2_n(DualSkewProduct(this->l23, DualNumber(phi_2.real(), alpha) + flip));
                    auto correction = this->l34.acos3(s.inverse() * this->l12, m2_n.inverse() * this->l12);

                    solutions.push_back({
                        phi_1 + phi_3 + DualNumber(phi_off, translation_off),
                        DualNumber(phi_2.real(), alpha) + flip,
                        correction
                    });
                }
                //generic
            } else {
                DualFrame m2(DualSkewProduct(this->l23, phi_2));
                auto phi_1 = (-this->l12).acos3(s * this->l34, m2 * this->l34);
                auto phi_3 = this->l34.acos3(s.inverse() * this->l12, m2.inverse() * this->l12);
                solutions.push_back({phi_1, phi_2, phi_3});
            }
        }
    }

    return solutions;
}
