//
// Created by sba on 06.07.21.
//

#include <iostream>
#include <numeric>

#include "dual_number.h"
#include "vector.h"
#include "unit_line.h"
#include "dual_frame.h"

#include "ccc.h"

#include "random.h"

using namespace DualNumberAlgebra;

DualNumber toDeg(const DualNumber &a) {
    return DualNumber(a.real()/M_PI*180.0, a.dual());
}

struct DataSet {
    UnitLine a;
    UnitLine b;
    UnitLine c;
    DualFrame zp;
    DualFrame goal;

    DataSet(UnitLine a, UnitLine b, UnitLine c, DualFrame zp, DualFrame goal)
    : a(std::move(a)), b(std::move(b)), c(std::move(c)), zp(std::move(zp)), goal(std::move(goal)) {}
};

// Could be done also as a config file but that would have to be read somehow
// I want to reduce the dependencies and this is only meant as examples especially for development

// Create an named example with three lines, a zero posture and a goal frame
std::unordered_map<std::string, DataSet> examples = {
        {"coincide", DataSet(
                UnitLine(
                        UnitDirectionVector(0.0, 0.0, 1.0),
                        PointVector(0.0, 0.0, 0.0)
                ),
                UnitLine(
                        UnitDirectionVector(1.0, 0.0, 0.0),
                        PointVector(0.0, -1.0, 0.0)
                ),
                UnitLine(
                        DirectionVector(0.0, 1.0, 0.0),
                        PointVector(10.0, 0.0, 8.0)
                ),
                DualFrame(
                        RotationMatrix(M_PI_2, M_PI_2, 0),
                        PointVector(10, 5, 7)
                ),
                DualFrame(
                        RotationMatrix(0, 0, 0),
                        PointVector(15, 0, 0)
                ))
        },
        {"coincide2", DataSet(
                UnitLine(
                        UnitDirectionVector(0.0, 0.0, 1.0),
                        PointVector(0.0, 0.0, 0.0)
                ),
                UnitLine(
                        UnitDirectionVector(1.0, 0.0, 0.0),
                        PointVector(0.0, -1.0, 0.0)
                ),
                UnitLine(
                        DirectionVector(0.0, 0.0, 1.0),
                        PointVector(5.0, 0.0, 0.0)
                ),
                DualFrame(
                        RotationMatrix(0, 0, 0),
                        PointVector(6, 1, 0)
                ),
                DualFrame(
                        RotationMatrix(0, 0, 0),
                        PointVector(1, 1, 0)
                ))
        },
        {"parallel", DataSet(
                UnitLine(
                        UnitDirectionVector(0.0, 0.0, 1.0),
                        PointVector(0.0, 0.0, 0.0)
                ),
                UnitLine(
                        UnitDirectionVector(1.0, 0.0, 0.0),
                        PointVector(0.0, -1.0, 0.0)
                ),
                UnitLine(
                        DirectionVector(0.0, -1.0, 1.0),
                        PointVector(1.0, 0.0, 0.0)
                ),
                DualFrame(
                        RotationMatrix(0,0,0),
                        PointVector(0,-1,4)
                ),
                DualFrame(
                        RotationMatrix(0,0,-M_PI_4),
                        PointVector(4,4,0)
                ))
        },
        {"simple-parallel", DataSet(
                UnitLine(
                        UnitDirectionVector(0.0, 0.0, 1.0),
                        PointVector(0.0, 0.0, 0.0)
                ),
                UnitLine(
                        UnitDirectionVector(1.0, 0.0, 0.0),
                        PointVector(0.0, -1.0, 0.0)
                ),
                UnitLine(
                        DirectionVector(0.0, 0, 1.0),
                        PointVector(2.0, 1.0, 0.0)
                ),
                DualFrame(
                        RotationMatrix(0,0,0),
                        PointVector(2,0,0)
                ),
                DualFrame(
                        RotationMatrix(0,0,0),
                        PointVector(4,4,4)
                ))
        },
        {"screwed", DataSet(
                UnitLine(
                        UnitDirectionVector(0.0, 0.0, 1.0),
                        PointVector(0.0, 0.0, 0.0)
                ),
                UnitLine(
                        UnitDirectionVector(1.0, 0.0, 0.0),
                        PointVector(0.0, -1.0, 0.0)
                ),
                UnitLine(
                        DirectionVector(0.0, 1.0, 0.0),
                        PointVector(1.0, 0.0, 1.0)
                ),
                DualFrame(
                        RotationMatrix(0,0,0),
                        PointVector(0,-1,4)
                ),
                DualFrame(
                        RotationMatrix(M_PI_2,0,-M_PI_4),
                        PointVector(4,4,0)
                ))
        },
};

void test_scenario(char * input) {
    try {
        // May throw an exception if the key is not existent
        DataSet d = examples.at(input);
        CCCMechanism ccc(d.a, d.b, d.c, d.zp);

        auto solutions = ccc.inverse(d.goal);

        std::cout << "To Pose: " << std::endl << d.goal << std::endl;

        for (auto solution : solutions) {
            std::cout << "Solution " << std::endl << "  " <<
                      solution.phi_1 << std::endl << "  " <<
                      solution.phi_2 << std::endl << "  " <<
                      solution.phi_3 << std::endl;
            std::cout << "yields to:" << std::endl << ccc.forward(solution) << std::endl;
        }
    } catch(const std::out_of_range &e) {
        throw std::logic_error("Example not found");
    }

}

int main(int argc, char **argv) {
    try {
        if (argc == 1) {
            throw std::logic_error("Need an example to run");
        }

        test_scenario(argv[1]);
    }
    catch(const std::logic_error &e) {
        std::cerr   << e.what()  << std::endl
                    << "Usage: " << std::endl
                    << argv[0] << " "
                    << std::accumulate(examples.cbegin()++, examples.cend(), examples.cbegin()->first, [] (auto a, auto b) -> std::string {
                        return a + " | " + b.first;
                    })
                    << std::endl;
    }
}


