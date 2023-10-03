cmake_minimum_required(VERSION 3.10)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE) # Look into this when have time.

# Require C++11 and disable compiler-specific extensions
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set the project name
project(transfer LANGUAGES C CXX)

if(APPLE)
  find_package(Python 3 COMPONENTS Interpreter Development NumPy REQUIRED)
  set(PYTHON_MAJOR $ENV{Python_VERSION_MAJOR})
  set(PYTHON_MINOR $ENV{Python_VERSION_MINOR})
  set(PYTHONLIBS_VERSION_STRING ${Python_VERSION})
  set(PYTHON_INCLUDE_DIR ${Python_INCLUDE_DIRS})
  set(PYTHON_LIBRARIES ${Python_LIBRARIES})
  message("Apple - Using Python:${Python_VERSION_MAJOR} - Libraries:${PYTHON_LIBRARIES} - IncludeDirs: ${PYTHON_INCLUDE_DIR}")
elseif(UNIX)
  find_package(PythonLibs 3 REQUIRED)
else()
  set(PYTHONLIBS_VERSION_STRING $ENV{PY_VERSION})
endif()

include_directories(${PYTHON_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIR})
# link numpy
link_libraries(Python::NumPy ${Python_LIBRARIES})

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
INCLUDE_DIRECTORIES(include src)

# populate matplotlib repository
include(FetchContent)
FetchContent_Declare(
    matplotlib
    GIT_REPOSITORY https://github.com/lava/matplotlib-cpp.git
    GIT_TAG        f23347fca25219d1c42cbb91608b5556814bf572
)
FetchContent_GetProperties(matplotlib)
if(NOT matplotlib_POPULATED)
    FetchContent_Populate(matplotlib)
endif()
include_directories(SYSTEM ${matplotlib_SOURCE_DIR})

add_subdirectory(src)

# Enable testing
# https://cmake.org/cmake/help/latest/command/add_test.html
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(CTest)
endif()

# Add test folder
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()

