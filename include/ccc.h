//
// Created by sba on 06.07.21.
//

#ifndef DAK_ADJOINT_TRIGONOMETRY_H
#define DAK_ADJOINT_TRIGONOMETRY_H

#include <vector>

#include "dual_number.h"
#include "types/matrix6.h"
#include "types/pluecker.h"

struct Configuration {
    DualNumberAlgebra::DualNumber phi_1;
    DualNumberAlgebra::DualNumber phi_2;
    DualNumberAlgebra::DualNumber phi_3;
};

struct CCCMechanism {
    Pluecker l12;
    Pluecker l23;
    Pluecker l34;
    AdjungateMatrix zero_posture;

    CCCMechanism(const Pluecker &l12, const Pluecker &l23, const Pluecker &l34, const AdjungateMatrix &zero_posture) noexcept;

    AdjungateMatrix forward(const Configuration &config) const noexcept;
    std::vector<Configuration> inverse(const AdjungateMatrix &pose) const;
};

#endif //DAK_ADJOINT_TRIGONOMETRY_H
