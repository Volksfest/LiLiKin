cmake_minimum_required(VERSION 3.10)

project("Dual Algebra Kinematics" VERSION 0.1)

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

# TODO check if googletest is pulled

add_subdirectory(googletest)
add_executable(do_tests
        src/dual_number.cpp
        test/dual_number_test.cpp)
target_include_directories(do_tests PRIVATE include)
target_link_libraries(do_tests gtest_main)

### Docs

# TODO only create docs, if Doxygen is found!

find_package(Doxygen)

set(DOXYGEN_PROJECT_LOGO "icon.svg")
set(DOXYGEN_EXTRA_PACKAGES "amsmath")
set(DOXYGEN_GENERATE_TODOLIST "YES")
doxygen_add_docs(docs
        ${PROJECT_SOURCE_DIR}/include
        COMMENT "Generate HTML Docs..."
        )
