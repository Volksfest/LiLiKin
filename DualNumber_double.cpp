//
// Created by sba on 17.09.20.
//

// Only needed because trigonometric only makes sense for doubles

#include "dual_number.h"

DualNumberAlgebra::DualNumber_template<double> DualNumberAlgebra::sqrt(const DualNumberAlgebra::DualNumber_template<double> & phi) noexcept{
    const double &sr = std::sqrt(phi.getReal());
    const double &d = phi.getDual();
    return DualNumberAlgebra::DualNumber_template<double>(sr, 0.5 * d / sr);
}

DualNumberAlgebra::DualNumber_template<double> DualNumberAlgebra::sin(const DualNumberAlgebra::DualNumber_template<double> & phi) noexcept{
    const double &r = phi.getReal();
    const double &d = phi.getDual();
    return DualNumberAlgebra::DualNumber_template<double>(std::sin(r), d * std::cos(r));
}

DualNumberAlgebra::DualNumber_template<double> DualNumberAlgebra::asin(const DualNumberAlgebra::DualNumber_template<double> & w) noexcept {
    const double &r = w.getReal();
    const double nr = std::asin(r);
    const double &d = w.getDual();
    return DualNumberAlgebra::DualNumber_template<double>(nr, d / std::cos(nr));
}

DualNumberAlgebra::DualNumber_template<double> DualNumberAlgebra::cos(const DualNumberAlgebra::DualNumber_template<double> & phi) noexcept {
    const double &r = phi.getReal();
    const double &d = phi.getDual();
    return DualNumberAlgebra::DualNumber_template<double>(std::cos(r), -d * std::sin(r));
}

DualNumberAlgebra::DualNumber_template<double> DualNumberAlgebra::acos(const DualNumberAlgebra::DualNumber_template<double> & w) noexcept {

    const double &r = w.getReal();
    const double nr = std::acos(r);
    const double &d = w.getDual();
    //return DualNumberAlgebra::DualNumber<double>(nr, -d / std::sin(nr));
    return DualNumberAlgebra::DualNumber_template<double>(nr, -d / std::sqrt(1-r*r));
}

DualNumberAlgebra::DualNumber_template<double> DualNumberAlgebra::tan(const DualNumberAlgebra::DualNumber_template<double> & phi) noexcept {
    const double &r = phi.getReal();
    const double &d = phi.getDual();
    const double cr = std::cos(r);
    return DualNumberAlgebra::DualNumber_template<double>(std::tan(r), d / cr / cr);
}

DualNumberAlgebra::DualNumber_template<double> DualNumberAlgebra::atan(const DualNumberAlgebra::DualNumber_template<double> & w) noexcept {
    const double &r = w.getReal();
    const double nr = std::atan(r);
    const double &d = w.getDual();
    const double cnr = std::cos(nr);
    return DualNumberAlgebra::DualNumber_template<double>(nr, -d * cnr * cnr);
}

DualNumberAlgebra::DualNumber_template<double> DualNumberAlgebra::atan2(const DualNumberAlgebra::DualNumber_template<double> &y, const DualNumberAlgebra::DualNumber_template<double> &x) noexcept {
    const double &xr = x.getReal();
    const double &yr = y.getReal();
    const double &xd = x.getDual();
    const double &yd = y.getDual();
    return DualNumberAlgebra::DualNumber_template<double>(std::atan2(yr, xr), (xr*yd - xd * yr) / ( xr * xr + yr * yr));
}

DualNumberAlgebra::DualNumber_template<double> DualNumberAlgebra::operator "" _s(long double dual) noexcept{
    return DualNumberAlgebra::DualNumber_template<double>(0.0, dual);
}

// Convenience for real part for dual angles
double DualNumberAlgebra::operator "" _d(long double degree) noexcept{
    return degree * M_PI / 180.0;
}
