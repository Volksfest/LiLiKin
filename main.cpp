//
// Created by sba on 29.06.20.
//

#include <iostream>
#include <eigen3/Eigen/Eigen>

#include "DualNumber.h"
#include "RobotAlgebra.h"

using namespace DualNumberAlgebra;
using namespace RobotAlgebra;

double r2d(double r) {
    return r/M_PI * 180.0;
}

int main() {

    auto ab = Pluecker<double>::fromPoints(Vec3<double>(1.0, 1.0, 0.0), Vec3<double>(2.0, 1.0, 0.0));

    auto phi = M_PI + 0.0_s;

    auto d = ab.REGG(phi);

    auto adj = ab.adjungatePluecker();
    std::cout << adj <<std::endl << std::endl;

    auto ws = waschmaschine(adj);
    std::cout << ws <<std::endl << std::endl;

    auto sw = sandwich(ws);
    std::cout << sw <<std::endl << std::endl;

    std::cout << "sin: " << std::endl << adjungateDualAngle(sin(phi)) <<std::endl << std::endl;
    std::cout << "cos: " << std::endl << adjungateDualAngle(cos(phi))  << std::endl << std::endl;

    std::cout << ws + sw << std::endl << std::endl;


    Mat3<double> rot;
    Vec3<double> trans;
    decompose(d, rot, trans);
    std::cout << d << std::endl << std::endl;

    std::cout << rot << std::endl <<std::endl<< trans << std::endl;

    return 0;
}
