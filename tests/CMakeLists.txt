cmake_minimum_required(VERSION 3.10)

include_directories(${TEST_SRC_DIR}/src
                    ${TEST_SRC_DIR}/include)

# TESTING
add_executable(TestLapack EXCLUDE_FROM_ALL test_lapack.cpp)
add_test(NAME TestLapack COMMAND TestLapack)

# Create a new target to build and run tests
# This ensure make all doesn't make every test at each compilation stage
add_custom_target(build_tests)
add_custom_target(run_tests COMMAND ${CMAKE_CTEST_COMMAND})
add_dependencies(run_tests build_tests)
add_dependencies(build_tests TestLapack)

