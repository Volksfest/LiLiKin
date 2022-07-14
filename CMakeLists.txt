cmake_minimum_required(VERSION 3.10)

project("Library for Line Kinematics" VERSION 0.2)

set(CMAKE_CXX_STANDARD 17)

### Project

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

set(SOURCE
        src/ccc.cpp

        src/base/dual_number.cpp
        src/base/vector.cpp
        src/base/matrix3.cpp

        src/screws/screw.cpp
        src/screws/screw_cos3.cpp
        src/screws/unit_line.cpp

        src/embedded_types/dual_embedded_matrix.cpp
        src/embedded_types/dual_frame.cpp
        src/embedded_types/dual_skew.cpp
        src/embedded_types/dual_skew_product.cpp

        src/util/random.cpp
        )

add_library(lilikin
        ${SOURCE}
        )
set_target_properties(lilikin PROPERTIES PUBLIC_HEADER
        "include/ccc.h;include/dual_number.h;include/vector.h;include/matrix3.h;include/screw.h;include/unit_line.h;include/dual_embedded_matrix.h;include/dual_frame.h;include/dual_skew.h;include/dual_skew_product.h;include/random.h;include/precision.h;include/lilikin.h")
target_include_directories(lilikin PRIVATE include)
target_link_libraries(lilikin Eigen3::Eigen)

### Example

add_executable(exlikin
    example.cpp)
target_include_directories(exlikin PRIVATE include)
target_link_libraries(exlikin PRIVATE lilikin)

### Shared Lib (optional)

# Needed so the install target only installs the shared lib if it is build
SET(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY TRUE)

# Could use the same precompiled objects with the static but I like the static library to be non PIC
add_library(lilikin_shared SHARED EXCLUDE_FROM_ALL
        ${SOURCE})
target_include_directories(lilikin_shared PRIVATE include)
target_link_libraries(lilikin_shared Eigen3::Eigen)

### Install

install(TARGETS lilikin lilikin_shared exlikin
        ARCHIVE
        DESTINATION lib
        LIBRARY
        DESTINATION lib
        OPTIONAL
        RUNTIME
        DESTINATION bin
        PUBLIC_HEADER
        DESTINATION include/lilikin
        )

### Test

if(EXISTS "${PROJECT_SOURCE_DIR}/googletest/CMakeLists.txt")
    add_subdirectory(googletest EXCLUDE_FROM_ALL)
    add_executable(lilikin_tests
            test/dual_number_test.cpp
            test/screws_test.cpp
            test/ccc_test.cpp
            test/random_test.cpp
            )
    target_include_directories(lilikin_tests PRIVATE include)
    target_link_libraries(lilikin_tests gtest_main lilikin)

    message(STATUS "Googletest cloned. Will also create tests. After build, run them with:")
    message(STATUS "    ./lilikin_tests")
else()
    message(STATUS "Googletest not cloned. Tests will be omitted!")
endif()

### Docs

find_package(Doxygen QUIET)

if(${DOXYGEN_FOUND})
    set(DOXYGEN_PROJECT_LOGO "doc/icon.svg")
    set(DOXYGEN_EXTRA_PACKAGES amsmath amssymb)
    set(DOXYGEN_GENERATE_TODOLIST "YES")
    set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "README.md")
    set(DOXYGEN_IMAGE_PATH ${PROJECT_SOURCE_DIR})
    set(DOXYGEN_SHOW_FILES "NO")
    doxygen_add_docs(docs
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/README.md
        COMMENT "Generate HTML Docs..."
        )
    message(STATUS "Doxygen found. Docs can be created with:")
    message(STATUS "    make docs")
elseif()
    message(STATUS "Doxygen not found. Docs cannot be build!")
endif()

### Additional messages

message(STATUS "If you need a shared object you will have to build it separately:")
message(STATUS "    make lilikin_shared")