cmake_minimum_required(VERSION 3.15)

project(cpy_credit LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O3 -Wall")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# --- Source Files ---
set(SOURCES
    src/credit/cpy_credit.cpp
    src/credit/risk_neutral.cpp
    src/credit/stvol.cpp
)

# Explicit inclusions for Apple Silicon
set(CMAKE_PREFIX_PATH "/opt/homebrew/")
include_directories("/opt/homebrew/include/")

# --- Find Dependencies ---
set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 REQUIRED)
pybind11_add_module(cpy_credit ${SOURCES})

find_package(NLopt REQUIRED)

# Find Boost Test framework, essential for CTest integration
find_package(Boost REQUIRED COMPONENTS unit_test_framework)


# Apply common definitions and linking for the library
target_compile_definitions(cpy_credit PRIVATE VERSION_INFO=0.1.0)
target_link_libraries(cpy_credit PRIVATE NLopt)


# Copy the built library to src/pyvol folder for Python packaging
add_custom_command(TARGET cpy_credit
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_FILE:cpy_credit> ${CMAKE_SOURCE_DIR}/src/pyvol/$<TARGET_FILE_NAME:cpy_credit>
    COMMENT "Copying library to pyvol source directory."
)

# --- Testing ---
enable_testing()

# Test Executable 1: Heston Test Suite (now integrated with Boost)
add_executable(HestonTest 
    tests/cpp_tests/HestonTest.cpp
    src/credit/stvol.cpp
    src/credit/risk_neutral.cpp
)
# Link test executable against NLopt and Boost::unit_test library
target_link_libraries(HestonTest PRIVATE NLopt Boost::unit_test_framework) 

# Register the test with CTest
add_test(HestonTest HestonTest)

# Test Executable 2: Secant Test
add_executable(SecantTest 
    tests/cpp_tests/SecantTest.cpp
)
target_link_libraries(SecantTest PRIVATE NLopt Boost::unit_test_framework) 
add_test(SecantTest SecantTest)
