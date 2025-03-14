cmake_minimum_required(VERSION 3.29)

project(cpy_credit LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O3 -Wall")

set(SOURCES
    src/credit/cpy_credit.cpp
    src/credit/risk_neutral.cpp
    src/credit/stvol.cpp
)

# Explicit inclusions for Apple Silicon
set(CMAKE_PREFIX_PATH "/opt/homebrew/")
include_directories("/opt/homebrew/include/")

find_package(pybind11 REQUIRED)
pybind11_add_module(cpy_credit ${SOURCES})

find_package(NLopt REQUIRED)

target_compile_definitions(cpy_credit PRIVATE VERSION_INFO=0.1.0)
target_link_libraries(cpy_credit PRIVATE NLopt)

# Copy the built library to src/pyvol folder for the purposes of python-packaging
add_custom_command(TARGET cpy_credit
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E rename
        $<TARGET_FILE:cpy_credit> ${CMAKE_SOURCE_DIR}/src/pyvol/$<TARGET_FILE_NAME:cpy_credit>
    COMMENT "Copying library to pyvol source directory."
)

# Testing
enable_testing()

add_executable(HestonTest 
    tests/cpp_tests/HestonTest.cpp
    src/credit/stvol.cpp
    src/credit/risk_neutral.cpp
)
target_link_libraries(HestonTest PRIVATE NLopt)
add_test(HestonTest HestonTest)

add_executable(SecantTest 
    tests/cpp_tests/SecantTest.cpp
)
add_test(SecantTest SecantTest)
