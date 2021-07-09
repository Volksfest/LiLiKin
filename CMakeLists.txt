cmake_minimum_required(VERSION 3.10)
project(dak)

set(CMAKE_CXX_STANDARD 17)

### Project

add_executable(dak
        src/main.cpp

        src/dual_number.cpp
        src/trigonometry_helper.cpp
        src/ccc.cpp

        src/types/pluecker.cpp
        src/types/vector.cpp
        src/types/matrix3.cpp
        src/types/matrix6.cpp
        src/types/homogenous.cpp

        )
target_include_directories(dak PRIVATE include)

### Test

add_subdirectory(googletest)

add_executable(test_dual_number
        src/dual_number.cpp
        test/dual_number_test.cpp)
target_include_directories(test_dual_number PRIVATE include)
target_link_libraries(test_dual_number gtest_main)
