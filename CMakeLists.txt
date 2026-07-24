cmake_minimum_required(VERSION 3.15)

project(cpy_credit LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(APPLE)
  execute_process(
    COMMAND brew --prefix
    OUTPUT_VARIABLE HOMEBREW_PREFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(HOMEBREW_PREFIX)
    list(APPEND CMAKE_PREFIX_PATH "${HOMEBREW_PREFIX}")
    list(APPEND CMAKE_INCLUDE_PATH "${HOMEBREW_PREFIX}/include")
    list(APPEND CMAKE_LIBRARY_PATH "${HOMEBREW_PREFIX}/lib")
  endif()
endif()

find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)

# Locate nanobind's CMake config shipped inside the installed Python package.
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE nanobind_ROOT
)
find_package(nanobind CONFIG REQUIRED)

find_package(NLopt REQUIRED)
find_package(Boost REQUIRED COMPONENTS unit_test_framework)

nanobind_add_module(cpy_credit
    src/credit/cpy_credit.cpp
    src/credit/risk_neutral.cpp
    src/credit/stvol.cpp
)

target_compile_definitions(cpy_credit PRIVATE VERSION_INFO=0.1.0)

target_include_directories(cpy_credit PRIVATE ${NLOPT_INCLUDE_DIRS})
if(TARGET NLopt::nlopt)
  target_link_libraries(cpy_credit PRIVATE NLopt::nlopt)
elseif(TARGET NLOPT::nlopt)
  target_link_libraries(cpy_credit PRIVATE NLOPT::nlopt)
else()
  target_link_libraries(cpy_credit PRIVATE ${NLOPT_LIBRARIES})
endif()

add_custom_command(TARGET cpy_credit
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_FILE:cpy_credit> ${CMAKE_SOURCE_DIR}/src/pyvol/$<TARGET_FILE_NAME:cpy_credit>
    COMMENT "Copying library to pyvol source directory."
)

enable_testing()

add_executable(HestonTest
    tests/cpp_tests/HestonTest.cpp
    src/credit/stvol.cpp
    src/credit/risk_neutral.cpp
)

target_include_directories(HestonTest PRIVATE ${NLOPT_INCLUDE_DIRS})
target_link_libraries(HestonTest PRIVATE Boost::unit_test_framework)
if(TARGET NLopt::nlopt)
  target_link_libraries(HestonTest PRIVATE NLopt::nlopt)
elseif(TARGET NLOPT::nlopt)
  target_link_libraries(HestonTest PRIVATE NLOPT::nlopt)
else()
  target_link_libraries(HestonTest PRIVATE ${NLOPT_LIBRARIES})
endif()
add_test(NAME HestonTest COMMAND HestonTest)

add_executable(RootSolverTest
    tests/cpp_tests/RootSolverTest.cpp
)

target_include_directories(RootSolverTest PRIVATE ${NLOPT_INCLUDE_DIRS})
target_link_libraries(RootSolverTest PRIVATE Boost::unit_test_framework)
if(TARGET NLopt::nlopt)
  target_link_libraries(RootSolverTest PRIVATE NLopt::nlopt)
elseif(TARGET NLOPT::nlopt)
  target_link_libraries(RootSolverTest PRIVATE NLOPT::nlopt)
else()
  target_link_libraries(RootSolverTest PRIVATE ${NLOPT_LIBRARIES})
endif()
add_test(NAME RootSolverTest COMMAND RootSolverTest)
