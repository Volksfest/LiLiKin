cmake_minimum_required(VERSION 3.10)
project(dak)

set(CMAKE_CXX_STANDARD 17)

add_executable(dak
        DualNumber_double.cpp
        RobotAlgebraTypes.cpp
        math.cpp
        main.cpp
        )
