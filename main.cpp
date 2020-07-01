//
// Created by sba on 29.06.20.
//

#include <iostream>
#include <eigen3/Eigen/Eigen>

#include "DualNumber.h"
#include "RobotAlgebra.h"

using namespace DualNumberAlgebra;
using namespace RobotAlgebra;

template<class T>
T filter_zero(const T &a, double epsilon = 1e-8) {
    return (a.array().abs() > epsilon).select(a, 0.0);
}

int main() {

    Pluecker<double> ab = create_pluecker_from_points(
            Vec3<double>(1.0, 1.0, 0.0),
            Vec3<double>(2.0, 1.0, 0.0));

    DualNumber<double> phi = 180.0_d + 2.0_s;

    Mat6<double> d = REGG(ab, phi);
    Mat4<double> hom = convert_mat6_to_mat4(d);

    std::cout << "Adjungierte Displacement: " << std::endl << filter_zero(d) << std::endl << std::endl;
    std::cout << "Rot:" << std::endl << filter_zero(get_rot(hom)) << std::endl <<std::endl<< "Trans:" << std::endl << filter_zero(get_trans(hom)) << std::endl;
    std::cout << "Hom:" << std::endl << filter_zero(hom) << std::endl << std::endl;
    std::cout << "Adjungierte Displacement: " << std::endl << filter_zero(convert_mat4_to_mat6(hom)) << std::endl << std::endl;

    return 0;
}
