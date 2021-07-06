//
// Created by sba on 29.06.20.
//

#include <iostream>

#include "DualNumber.h"
//#include "RobotAlgebra.h"
#include "RobotAlgebraTypes.h"

#include "math.h"

using namespace RobotAlgebra;
using DualNumberAlgebra::operator""_s;

template<class T>
T filter_zero(const T &a, double epsilon = 1e-8) {
    return (a.array().abs() > epsilon).select(a, 0.0);
}

DualNumber toDeg(const DualNumber &a) {
    return DualNumber(a.getReal()/M_PI*180.0, a.getDual());
}


int main() {
    Pluecker a = Pluecker::fromDirection(
            Vector(0.0, 0.0, 1.0),
            Vector(0.0, 0.0, 0.0)
    );

    Pluecker b = Pluecker::fromDirection(
            Vector(1.0, 0.0, 0.0),
            Vector(0.0, 0.0, 0.0)
    );

    Pluecker c = Pluecker::fromDirection(
            Vector(0.0, 0.0, 1.0),
            Vector(0.0, 0.0, 0.0)
    );

    auto from_hom = HomogenousMatrix::from(
            RotationMatrix(0,0,0),
            Vector(2,0,4)
    );

    auto from(AdjungateMatrix::from(from_hom) );

    auto to_hom = HomogenousMatrix::from(
            RotationMatrix(0,0,0),
            Vector(-5,5,8)
    );

    auto to = AdjungateMatrix::from(to_hom);

    DualNumber phi1, phi2, phi3;
    solve(to, from, a, b, c, phi1, phi2, phi3, false);

    std::cout <<std::endl <<
            toDeg(phi1) << std::endl <<
            toDeg(phi2) << std::endl <<
            toDeg(phi3) << std::endl << std::endl;


    HomogenousMatrix res(HomogenousMatrix::from(
            AdjungateMatrix(a.transform(phi1).data * b.transform(phi2).data * c.transform(phi3).data)
    ));

    std::cout << res.data * from_hom.data << std::endl << std::endl;
    std::cout << filter_zero(to_hom.data) << std::endl << std::endl;

}
