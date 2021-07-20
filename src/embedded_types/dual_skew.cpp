//
// Created by sba on 19.07.21.
//

#include "embedded_types/dual_skew.h"
#include "screws/screw.h"
#include "base/vector.h"

DualSkew::DualSkew(const SkewMatrix &real, const SkewMatrix &dual) noexcept: DualEmbeddedMatrix(real, dual) {}

DualSkew::DualSkew(const UnitLine &line) noexcept: DualSkew(SkewMatrix(line.n()), SkewMatrix(line.m())) {}
