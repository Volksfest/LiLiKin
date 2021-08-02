//
// Created by sba on 06.07.21.
//

#ifndef DAK_ADJOINT_TRIGONOMETRY_H
#define DAK_ADJOINT_TRIGONOMETRY_H

#include <vector>
#include <tuple>

#include "base/dual_number.h"
#include "screws/unit_line.h"
#include "embedded_types/dual_frame.h"

/**
 * \brief The configuration of an CCC mechanism
 *
 * The configuration is given by the values of each joint.
 */
struct Configuration {
    DualNumberAlgebra::DualNumber phi_1; //!< Value of the first joint
    DualNumberAlgebra::DualNumber phi_2; //!< Value of the second joint
    DualNumberAlgebra::DualNumber phi_3; //!< Value of the third joint
};

/**
 * \brief A CCC mechanism defined by three lines and a zero posture frame pointing the endeffector with zero valued joints
 *
 */
struct CCCMechanism {
    UnitLine l12; //!< Line for the first C joint
    UnitLine l23; //!< Line for the second C joint
    UnitLine l34; //!< Line for the third C joint
    DualFrame zero_posture; //!< Zero posture frame in 6x6 matrix representation

    /**
     * \brief Simple constructor for the CCC mechanism
     * @param l12 First C joint
     * @param l23 Second C joint
     * @param l34 Third C joint
     * @param zero_posture Endeffector pose in zeroed joint values
     */
    CCCMechanism(const UnitLine &l12, const UnitLine &l23, const UnitLine &l34, const DualFrame &zero_posture) noexcept;

    /**
     * \brief Forward kinematics with PoE
     * @param config The joint configuration to calculate the endeffector pose
     * @return The endeffector pose
     */
    DualFrame forward(const Configuration &config) const noexcept;

    /**
     * \brief Verbose Forward kinematics with PoE
     *
     * This one also returns the intermediate line position yielding from the transformations of the lines before.
     * The first line is omitted as its parameter will not change.
     *
     * @param config The joint configuration to calculate the endeffector pose
     * @return A tuple containing endeffector pose an the intermediate lines from the FK
     */
    std::tuple<DualFrame, UnitLine, UnitLine> forward_verbose(const Configuration &config) const;

    /**
     * \brief The inverse kinematics
     *
     * The solution is not unique thus a list of Configuration as a solution
     * \except std::logic_error If no solution is possible
     * @param pose The frame to reach
     * @return A list with possible configurations
     */
    std::vector<Configuration> inverse(const DualFrame &pose) const;
};

#endif //DAK_ADJOINT_TRIGONOMETRY_H
