cmake_minimum_required(VERSION 3.10)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE) # Look into this when have time.

# Require C++11 and disable compiler-specific extensions
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set the project name
project(rad_transfer LANGUAGES C CXX)

# Prevent in source builds
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/CMake)
include(PreventInSourceBuilds)

# find and link eigen
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
link_libraries(Eigen3::Eigen)

# Load user settings before the defaults - this way the defaults will not
# overwrite the user set options. If the user has not set all options, we still
# have the defaults.
# set(USER_CONFIG "${CMAKE_CURRENT_SOURCE_DIR}/config/user.cmake" CACHE PATH
#   "Path to optional user configuration file.")
# message(STATUS "(optional) USER_CONFIG = ${USER_CONFIG}")
# include("${USER_CONFIG}" OPTIONAL)
# include("${CMAKE_CURRENT_SOURCE_DIR}/config/defaults.cmake")

# Allow overwriting of the compiler by setting CXX/MPICXX on the command line or
# in user.cmake.
if (NOT CMAKE_CXX_COMPILER)
  if (CXX)
    set(CMAKE_CXX_COMPILER ${CXX})
    # Avoid some issues when CXX is defined
    unset(CXX)
    unset(CXX CACHE)
  endif()
  if (MPICXX)
    # In parallel MPICXX takes precedence, if defined.
    set(CMAKE_CXX_COMPILER ${MPICXX})
    # Setting the variables below circumvents autodetection, see FindMPI.cmake.
    set(MPI_CXX_INCLUDE_PATH "")
    set(MPI_CXX_LIBRARIES "")
  endif()
endif()

### Include Directories
include_directories(include)
add_subdirectory(src)

# Enable testing
# https://cmake.org/cmake/help/latest/command/add_test.html
include(CTest)
enable_testing()
add_subdirectory(tests)