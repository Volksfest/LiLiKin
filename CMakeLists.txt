cmake_minimum_required(VERSION 3.10)

project("Dual Algebra Kinematics" VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)

### Project

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_library(dak
        src/ccc.cpp

        src/base/dual_number.cpp
        src/base/vector.cpp
        src/base/matrix3.cpp

        src/screws/screw.cpp
        src/screws/screw_cos3.cpp
        src/screws/unit_screw.cpp
        src/screws/line.cpp
        src/screws/unit_line.cpp

        src/embedded_types/dual_embedded_matrix.cpp
        src/embedded_types/dual_frame.cpp
        src/embedded_types/dual_skew.cpp
        src/embedded_types/dual_skew_product.cpp
        )
target_include_directories(dak PRIVATE include)
target_link_libraries(dak Eigen3::Eigen)

### Example
add_executable(dak-example
    applications/example.cpp)
target_include_directories(dak-example PRIVATE include)
target_link_libraries(dak-example PRIVATE dak)


### Test

if(EXISTS "${PROJECT_SOURCE_DIR}/googletest/CMakeLists.txt")
    add_subdirectory(googletest)
    add_executable(do_tests
            test/dual_number_test.cpp
            test/screws_test.cpp
            test/ccc_test.cpp
            )
    target_include_directories(do_tests PRIVATE include)
    target_link_libraries(do_tests gtest_main dak)

    message(STATUS "Googletest cloned. Will also create tests")
else()
    message(NOTICE "Googletest not cloned. Tests will be omitted!")
endif()

### Docs

find_package(Doxygen QUIET)

if(${DOXYGEN_FOUND})
    set(DOXYGEN_PROJECT_LOGO "icon.svg")
    set(DOXYGEN_EXTRA_PACKAGES "amsmath")
    set(DOXYGEN_GENERATE_TODOLIST "YES")
    doxygen_add_docs(docs
        ${PROJECT_SOURCE_DIR}/include
        COMMENT "Generate HTML Docs..."
        )
    message(STATUS "Doxygen found. Docs can be created with:")
    message(STATUS "    make docs")
elseif()
    message(NOTICE "Doxygen not found. Docs cannot be build!")
endif()
