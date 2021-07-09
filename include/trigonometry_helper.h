//
// Created by sba on 06.07.21.
//

#ifndef DAK_TRIGONOMETRY_HELPER_H
#define DAK_TRIGONOMETRY_HELPER_H

#include "types.h"

using namespace DualNumberAlgebra;

DualNumber acos3(const Pluecker &a, const Pluecker &b, const Pluecker &n) noexcept;
std::vector<DualNumber> solve_trig(const DualNumber &a,const DualNumber &b, const DualNumber &c);

#endif //DAK_TRIGONOMETRY_HELPER_H
