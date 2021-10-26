//
// Created by sba on 03.09.21.
//

#ifndef DUAL_ALGEBRA_KINEMATICS_PRECISION_H
#define DUAL_ALGEBRA_KINEMATICS_PRECISION_H

#include <stdlib.h>

class Compare {
private:
    double epsilon = 1e-12;

public:
    inline double get_precision() { return this->epsilon; }
    inline void set_precision(double precision) { if (precision > 0) this->epsilon = precision; }

    inline bool equal(double a, double b) { return (abs(a-b) < this->epsilon); }

// Add a lot of syntactic sugar
    inline static bool is_equal(double a, double b) {
        return Compare::instance().equal(a,b);
    }

    inline static bool is_zero(double a) {
        return Compare::instance().equal(a, 0.0);
    }

// Make class a Singleton
private:
    Compare() = default;
public:
    static Compare & instance(){
        static Compare instance;
        return instance;
    }

    Compare(Compare const &) = delete;
    void operator=(Compare const &) = delete;
};

#endif //DUAL_ALGEBRA_KINEMATICS_PRECISION_H
