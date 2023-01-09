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
    // Reformulate the pose with the zero posture such that
    // S = M1 * M2 * M3
    // instead of
    // pose = M1 * M2 * M3 * zero_posture
    auto s = pose * this->zero_posture.inverse();

    // Calculate the projections used for the generalized rodriguez formula
    DualSkew crossterm(this->l23);
    DualEmbeddedMatrix uniterm = - crossterm * crossterm;
    DualEmbeddedMatrix squareterm = DualEmbeddedMatrix(1) - uniterm;

    // Calculate parameters regarding the rodriguez formula
    DualNumber a = this->l12 * (uniterm * this->l34);
    DualNumber b = this->l12 * (crossterm * this->l34);
    DualNumber c = this->l12 * ((s - squareterm) * this->l34); // TODO simplify

    // Calculate phi_2 as the trigonometric solutions of a cos + b sin = c
    std::vector<DualNumber> phi2_solutions = solve_trigonometric_equation(a, b, c);

    // Check the line relation between Line 1 and the final Line 3
    // This will result in annoying special cases
    LineRelation parallelity = this->l12.get_relation_to(s * this->l34);

    // Premodifier
    // As we are manipulating phi2_solutions, we use classical array access instead of nice fancy for-ranges :(
    // Maybe not the best way... TODO improve?
    for (int i = 0; i < phi2_solutions.size(); i++) {
        // 180° rotations have to be considered separately
        if (parallelity == LineRelation::ANTI_COINCIDE || parallelity== LineRelation::ANTI_PARALLEL) {
            phi2_solutions[i] +=  M_PI;
        }

        // atan2 with a solution of primal zero is not solveable right now
        // This happens in parallel cases and thus needs extra care
        if (parallelity == LineRelation::PARALLEL || parallelity== LineRelation::ANTI_PARALLEL) {
            DualFrame pre_rot(DualSkewProduct(this->l23, phi2_solutions[i].real()));
            auto l3i = pre_rot * this->l34;
            auto l3f = s * this->l34;

            double tri_b = abs(this->l12.get_distance(l3i).dual());
            double tri_c = abs(this->l12.get_distance(l3f).dual());

            // Here can actually be seen that to the first atan2 a second pure dual atan2 will be added/substracted
            // The "length" has to be the  dual part of the atan2 solution which cannot be retrieved.
            double p_2 = phi2_solutions[i].dual();
            double length = sqrt(p_2 * p_2 - tri_b * tri_b + tri_c * tri_c);

            // Make the actual two solutions of the only one in the parallel case
            phi2_solutions.push_back(phi2_solutions[i] + DualNumber(0, length));
            phi2_solutions[i] += DualNumber(0, -length);
            // Probably needed as the size increases. Not the best solution. It should really be improved...
            break;
        }
    }

    // Container for the solutions
    std::vector<Configuration> solutions;
    for (auto phi_2 : phi2_solutions) {
        // M2 can be calculated already and is the same in all cases
        DualFrame m2(DualSkewProduct(this->l23, phi_2));

        // Generic case
        if (parallelity == LineRelation::SKEW || parallelity == LineRelation::INTERSECT) {
            // Calculate the angles as the missing transformation around a single line
            // See paper: "The adjoint trigonometric representation of displacements
            // and a closed-form solution to the IKP of general 3C chains", Bongardt, ZAMM, 2019
            auto phi_1 = this->l12.acos3(m2 * this->l34, s * this->l34);
            auto phi_3 = this->l34.acos3(s.inverse() * this->l12, m2.inverse() * this->l12);
            solutions.push_back({phi_1, phi_2, phi_3});

        } else {
            // Coincide
            if (parallelity == LineRelation::ANTI_COINCIDE || parallelity == LineRelation::COINCIDE) {
                // Get an orthogonal line to check orientation of the frame
                // The orthogonality ensures, that it is not the rotation axis
                DirectionVector u(cross(this->l12.n(), this->l23.n()));
                UnitLine orthogonal(
                        u.normal(),
                        PointVector(0,0,0)
                );

                // Calculate the angles
                auto phi_1 = this->l12.acos3(m2 * orthogonal, s * orthogonal);
                // phi_3 is totally redundant as it gives the transformation of the coinciding line
                auto phi_3 = DualNumber();
                solutions.push_back({phi_1, phi_2, phi_3});

            // Parallel
            } else {
                auto phi_1 = this->l12.acos3(m2 * this->l34, s * this->l34);
                auto phi_3 = this->l34.acos3(s.inverse() * this->l12, m2.inverse() * this->l12);

                DualFrame pseudo_pose = this->forward({phi_1, phi_2, phi_3});
                auto d = (pose.p() - pseudo_pose.p()) * this->l12.n();
                solutions.push_back({phi_1 + DualNumber(0,d), phi_2, phi_3});
            }
        }
    }

    return solutions;
}
