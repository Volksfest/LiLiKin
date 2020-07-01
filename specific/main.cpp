#include <iostream>
#include <eigen3/Eigen/Eigen>

#include "DualNumber.h"
#include "RobotAlgebra.h"

using namespace DualNumberAlgebra;
using namespace RobotAlgebra;

void playaround1() {
    DualNumber a = DualNumber(3.0,2.0);
    std::cout << a.getReal() << " "<< a.getDual() << std::endl;
    std::cout << a << std::endl << std::endl;

    Eigen::Matrix<DualNumber, 6,6> m = Eigen::Matrix<DualNumber, 6, 6>::Random(6,6);
    std::cout << m << std::endl << "*" << std::endl;

    Eigen::Matrix<DualNumber, 6,6> n(6,6);
    n.topLeftCorner(3, 3)     = Eigen::Matrix<DualNumber, 3,3>::Zero(3, 3);
    n.topRightCorner(3, 3)     = Eigen::Matrix<DualNumber, 3,3>::Identity(3, 3);
    n.bottomLeftCorner(3, 3)     = Eigen::Matrix<DualNumber, 3,3>::Identity(3, 3);
    n.bottomRightCorner(3, 3)     = Eigen::Matrix<DualNumber, 3,3>::Zero(3, 3);
    std::cout << n << std::endl << "=" << std::endl;

    std::cout << m * n << std::endl;
}

int main() {

    auto ab = Pluecker::fromPoints(Vector3(1.0, 0.0, 0.0), Vector3(2.0, 2.0, 0.0));
    auto cd = Pluecker::fromPoints(Vector3(-1.0, 2.0, 0.0), Vector3(-1.0, 2.0, 1.0));
    auto ef = Pluecker::fromPoints(Vector3(-3.0, 0.0, 0.0), Vector3(-3.0, 3.0, 0.0));

    std::cout << ab << std::endl;

    /*
    std::cout << "ab:" << std::endl << ab.getDualLine() << std::endl << std::endl;
    std::cout << "cd:" << std::endl << cd.getDualLine() << std::endl << std::endl;
    std::cout << "ef:" << std::endl << ef.getDualLine() << std::endl << std::endl;

    std::cout << "ab*cd:" << std::endl << ab.getDualLine() * cd.getDualLine() << std::endl << std::endl;*/


    return 0;
}
