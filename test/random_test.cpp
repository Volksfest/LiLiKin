//
// Created by sba on 22.11.21.
//

#include "random.h"
#include <gtest/gtest.h>

// Just check, if they don't throw an exception
TEST(Random, Typecheck) {
    try {
        auto l = Random::SampleLine();
    } catch(...) {
        FAIL();
    }

    try{
        auto d = Random::SampleDirection();
    } catch(...) {
        FAIL();
    }

    try{
        auto o = Random::SampleOrientation();
    } catch(...) {
        FAIL();
    }

    try{
        auto f = Random::SampleFrame();
    } catch(...) {
        FAIL();
    }
}