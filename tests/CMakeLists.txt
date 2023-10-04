cmake_minimum_required(VERSION 3.10)

include_directories(${CMAKE_SOURCE_DIR}/src
                    ${CMAKE_SOURCE_DIR}/include)

# TESTING
add_executable(RandomTest EXCLUDE_FROM_ALL random_test.cpp)
add_test(NAME RandomTest COMMAND RandomTest)

add_executable(TestPlanckian EXCLUDE_FROM_ALL test_planckian.cpp ${rt_SOURCES})
add_test(NAME TestPlanckian COMMAND TestPlanckian)

# Create a new target to build and run tests
# This ensure make all doesn't make every test at each compilation stage
add_custom_target(build_tests)
add_custom_target(run_tests COMMAND ${CMAKE_CTEST_COMMAND})
add_dependencies(run_tests build_tests)
add_dependencies(build_tests RandomTest TestPlanckian)
