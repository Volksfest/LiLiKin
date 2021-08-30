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
    DualNumber c = this->l12 * ((s - squareterm) * this->l34);

    // Calculate phi_2 as the trigonometric solutions of a cos + b sin = c
    std::vector<DualNumber> phi2_solutions = solve_trigonometric_equation(a, b, c);

    // Check the line relation between Line 1 and the final Line 3
    // This will result in annoying special cases
    Parallelity parallelity = this->l12.is_parallel(s * this->l34);

    // Container for the solutions
    std::vector<Configuration> solutions;
    for (auto phi_2 : phi2_solutions) {
        // 180° rotations have to be considered separately
        if (parallelity == Parallelity::ANTI_COINCIDE || parallelity== Parallelity::ANTI_PARALLEL) {
            phi_2 = phi_2 + M_PI;
        }
        // M2 can be calculated already and is the same in all cases
        DualFrame m2(DualSkewProduct(this->l23, phi_2));

        // Generic case
        if (parallelity == Parallelity::SKEW || parallelity == Parallelity::INTERSECT) {
            // Calculate the angles as the missing transformation around a single line
            // See paper; TODO correct refs
            auto phi_1 = this->l12.acos3(m2 * this->l34, s * this->l34);
            auto phi_3 = this->l34.acos3(s.inverse() * this->l12, m2.inverse() * this->l12);
            solutions.push_back({phi_1, phi_2, phi_3});

        } else {
            // Coincide
            if (parallelity == Parallelity::ANTI_COINCIDE || parallelity == Parallelity::COINCIDE) {
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
                // Calculate the angles as usual first
                // These calculations will give the correct angles
                //   but unfortunately not the correct translations as they are not unique
                auto phi_1 = this->l12.acos3(m2 * this->l34, s * this->l34);
                auto phi_3 = this->l34.acos3(s.inverse() * this->l12, m2.inverse() * this->l12);

                // Calculate the only rotation displacements
                // As Line 1 and Line 3 will be parallel the rotations will be summed up in Line 1
                // This yields to different translational displacement but has no effect on the orientation
                DualFrame m1_combined(DualSkewProduct(this->l12,(phi_1 + phi_3).real()));
                DualFrame m2_rot(DualSkewProduct(this->l23, phi_2.real()));

                // Calculate a new zero posture which already has the correct orientation
                // Thus, the translational displacement of this pose has to be finded
                auto pre_pose = m1_combined * m2_rot * this->zero_posture;

                // We are using some kind of a Paden-Kahan problem
                // The problem is projected to an arbitrary plane (we just use the canonical anchor) orthogonal to Line 1/3
                PointVector anchor = this->l12.get_canonical_anchor();
                // The target and the new zero posture pose is projected to the plane
                PointVector p_target = this->l12.orthogonal_plane_projection(anchor,pose.p());
                PointVector p_begin = this->l12.orthogonal_plane_projection(anchor, (pre_pose).p());

                // We also define a new origin, which is also the canonical anchor
                PointVector o_o = this->l12.get_canonical_anchor();
                // We will have a second rotational center, which is the intersection of Line 3 with the plane
                // The intersection is calculated by projecting the origin to Line 3
                // The intersection is the closest point to any point on the plane
                PointVector o_r = (m1_combined * m2_rot * this->l34).point_project(o_o);

                // By calculating p_t = Rot_o(-phi) * Trans_n(alpha) * Rot(phi) * p_b
                // we receive the needed translational displacements inside the plane
                // Thus we need phi, the angle to rotate around Line 1 and counter-rotate around Line 3
                // and we need alpha, the translation around Line 2 (given by n)
                // The counter-rotation on Line 3 is needed to hold the correct orientation
                // By reducing the above formula we get:
                // p_t = p_b - o + Rot(phi) (o + alpha * n)
                // This has to be solved for alpha and phi

                // This p is acually p_t - p_b making the formula above easier
                Vector p = p_target - p_begin;
                // The rotation center by the new origin o_o
                Vector o = o_r - o_o;
                // The direction of Line 2 after the pre orientational displacements
                Vector n = (m1_combined * this->l23).n();
                // The orthogonal of the plane
                Vector w = this->l12.n();

                // Reformulation above with p yields to:
                // p + o = Rot(phi) (o + alpha * n)
                // The rotation will conserve any vector length, thus:
                // || p + o || = || o + alpha * n ||
                // alpha needs to be calculated such that the length are the same
                // This yields to a quadratic equation
                double minus_p_2 = -(n * o);
                double q = - (p * ( p + 2 * o));
                double sq = sqrt(minus_p_2 * minus_p_2 - q);
                // yielding to two solutions!
                std::vector<double> alphas = {minus_p_2 - sq, minus_p_2 + sq};

                // For each alpha we can now solve a Paden-Kahan problem of the first kind to get phi
                // v is the target projected vector for the Paden-Kahan problem
                Vector v = p + o;
                for (double alpha : alphas) {
                    // u is the starting projected vector for the Paden-Kahan problem which is dependend on alpha
                    Vector u = o + alpha * n;
                    // The result of the Paden-Kahan problem
                    double phi_off = atan2(w * cross(u,v), u * v);

                    // The last translation is given by the orthogonal translation to the plane
                    // The only translation outside the plane is given by the rotation around Line 2
                    // considering that and projecting all the points to Line 1, the offset can be calculated
                    auto tmp_translation = this->l12.point_project(pose.p()) - this->l12.point_project((m2 * this->zero_posture).p());
                    // The projection does not take into account in which direction the Line (actually spear) looks
                    // Thus it has to be considered additionally
                    double translation_off = tmp_translation.norm() * (tmp_translation * this->l12.n() > 0 ? 1.0 : -1.0);

                    DualNumber corrected_phi_2(phi_2.real(), alpha);
                    // Actually, it should be enough to set phi_3 as the counter-rotation of Line 1 being -phi_off
                    // But unfortunately, Through the parallelity there are additionally offsets not considered
                    // To make it easy, just calculate the simple acos3 between what we have and what we should have
                    DualFrame m2_n(DualSkewProduct(this->l23, corrected_phi_2));
                    auto corrected_phi_3 = this->l34.acos3(s.inverse() * this->l12, m2_n.inverse() * this->l12);

                    solutions.push_back({
                        phi_1 + phi_3 + DualNumber(phi_off, translation_off),
                        corrected_phi_2,
                        corrected_phi_3
                    });
                }
            }
        }
    }

    return solutions;
}
