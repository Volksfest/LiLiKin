//
// Created by sba on 22.11.21.
//

#include "random.h"

#include <random>
#include <functional>

double RandomNumber() {
    static std::random_device rd;  // Will be used to obtain a seed for the random number engine
    static std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    static std::uniform_real_distribution<> dis(0.0, 1.0);

    return dis(gen);
}

Vector Random::SampleVector() {
    return Vector(
            RandomNumber(),
            RandomNumber(),
            RandomNumber()
            );
}

// https://www.bogotobogo.com/Algorithms/uniform_distribution_sphere.php
UnitDirectionVector Random::SampleDirection() {
    double theta = 2 * M_PI * RandomNumber();
    double phi = acos(2 * RandomNumber() - 1.0);

    return UnitDirectionVector(
            cos(theta) * sin(phi),
            sin(theta) * sin(phi),
            cos(phi));
}

// See: Fast Random Rotation Matrices - James Arvo - DOI: https://doi.org/10.1016/B978-0-08-050755-2.50034-8
RotationMatrix Random::SampleOrientation() {
    double x1 = RandomNumber();
    double x2 = RandomNumber();
    double x3 = RandomNumber();

    double PI2 = 2 * M_PI;
    double c1 = cos(PI2*x1);
    double s1 = sin(PI2*x1);
    double sqr3 = sqrt(x3);

    Eigen::Matrix<double, 3, 3> mat;
    Eigen::Matrix<double, 3, 1> vec;
    mat << c1, s1, 0, -s1, c1, 0, 0, 0,1;
    vec << cos(PI2 * x2) * sqr3, sin(PI2 * x2) * sqr3, sqrt(1-x3);

    Eigen::Matrix<double, 3, 3> rot = -(Eigen::Matrix<double, 3, 3>::Identity(3,3) - 2 * vec * vec.transpose()) * mat;
    return RotationMatrix::RotationFromEigen(rot);
}

PointVector Random::SamplePoint(double cube_size) {
    return PointVector(cube_size * RandomNumber() * Random::SampleVector());
}

UnitLine Random::SampleLine() {
    return UnitLine(
            SampleDirection(),
            SamplePoint(5));
}

DualFrame Random::SampleFrame() {
    PointVector p = PointVector(5 * RandomNumber() * Random::SampleVector());
    return DualFrame(Random::SampleOrientation(), p);
}

DualNumberAlgebra::DualNumber Random::SampleDualNumber() {
    return DualNumberAlgebra::DualNumber(RandomNumber(), RandomNumber());
}

UnitLine Random::SampleRelatedLine(const UnitLine &a, LineRelation relation) {
    UnitDirectionVector n = SampleDirection();
    PointVector anchor = SamplePoint(5);

    if (relation == LineRelation::INTERSECT || relation == LineRelation::ANTI_COINCIDE || relation == LineRelation::COINCIDE) {
        anchor = a.get_canonical_anchor();
    }

    if (relation == LineRelation::PARALLEL || relation == LineRelation::COINCIDE) {
        n = a.n().normal();
    }

    if (relation == LineRelation::ANTI_COINCIDE || relation == LineRelation::ANTI_PARALLEL) {
        n = UnitDirectionVector(a.n() * -1);
    }

    return UnitLine(n, anchor);
}


std::pair<UnitLine, UnitLine> Random::SampleLinePair(LineRelation relation) {
    UnitLine a = SampleLine();

    return std::make_pair(a, SampleRelatedLine(a, relation));
}

std::tuple<UnitLine, UnitLine, UnitLine>
Random::SampleLineTriplet(LineRelation relation_a_to_b, LineRelation relation_a_to_ref) {
    UnitLine a = SampleLine();
    UnitLine b = SampleRelatedLine(a, relation_a_to_b);
    UnitLine ref = SampleRelatedLine(a, relation_a_to_ref);

    return std::make_tuple(a,b, ref);
}

std::tuple<UnitLine, UnitLine, UnitLine>
Random::SamplePlanarLineTriplet(bool intersecting_a_b, bool intersecting_a_ref) {
    UnitLine a = SampleLine();
    UnitLine b = SampleRelatedLine(a, intersecting_a_b? LineRelation::INTERSECT : LineRelation::SKEW);

    double alpha = RandomNumber() - 0.5;
    double beta = RandomNumber() - 0.5;

    return std::make_tuple(
            a,
            b,
            // create the third line which has a planar direction
            UnitLine(
                    // planarity is given by a random composition of the other both directions
                    DirectionVector(alpha * a.n() + beta * b.n()).normal(),
                    // if a shall intersect with reference than the same anchor is used
                    // otherwise a random anchor will be sampled
                    intersecting_a_ref?
                        a.get_canonical_anchor():
                        SamplePoint(5)
                    )
            );
}
