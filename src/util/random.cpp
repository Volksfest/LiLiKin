//
// Created by sba on 22.11.21.
//

#include "util/random.h"

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

UnitLine Random::SampleLine() {
    return UnitLine(
            Random::SampleDirection(),
            PointVector(5 * RandomNumber() * Random::SampleVector()));
}

DualFrame Random::SampleFrame() {
    PointVector p = PointVector(5 * RandomNumber() * Random::SampleVector());
    return DualFrame(Random::SampleOrientation(), p);
}

DualNumberAlgebra::DualNumber Random::SampleDualNumber() {
    return DualNumberAlgebra::DualNumber(RandomNumber(), RandomNumber());
}