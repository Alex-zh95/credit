cmake_minimum_required(VERSION 3.29)

project(cpy_credit LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O3 -Wall")

set(SOURCES
    src/credit/cpy_credit.cpp
    src/credit/risk_neutral.cpp
    src/credit/stvol.cpp
)

find_package(pybind11 REQUIRED)
pybind11_add_module(cpy_credit ${SOURCES})

# Additional target for NLopt library
find_package(NLopt REQUIRED)

target_compile_definitions(cpy_credit PRIVATE VERSION_INFO=0.1.0)

# Copy the built library to src/pyvol folder for the purposes of python-packaging
add_custom_command(TARGET cpy_credit
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E rename
        $<TARGET_FILE:cpy_credit> ${CMAKE_SOURCE_DIR}/src/pyvol/$<TARGET_FILE_NAME:cpy_credit>
    COMMENT "Copying library to pyvol source directory."
)

# Testing
enable_testing()

add_executable(HestonTest tests/cpp_tests/HestonTest.cpp src/credit/stvol.cpp)
target_link_libraries(HestonTest PRIVATE nlopt)
add_test(HestonTest HestonTest)

add_executable(SecantTest tests/cpp_tests/secant_test.cpp)
add_test(SecantTest SecantTest)
