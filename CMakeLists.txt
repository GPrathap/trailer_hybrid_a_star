cmake_minimum_required (VERSION 3.5)
project (lib_cpp_hybrid_a_star)

set(CMAKE_BUILD_TYPE "DEBUG")
set(CMAKE_CXX_FLAGS "-std=c++17")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g --debug-cpp --trace-expand")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_MATH_DISABLE_FLOAT128 -O1 -Wall -g --debug")

find_package(Eigen3 REQUIRED)
# find python libraries
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
find_package(PythonLibs 3.0 REQUIRED)
include_directories(${PYTHON3_INCLUDE_DIRS} ${NumPy_INCLUDE_DIRS})

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

include_directories( 
    SYSTEM 
    include
    ${Eigen3_INCLUDE_DIRS}
    ${matplotlib_SOURCE_DIR}
)

add_library(lib_cpp_hybrid_a_star src/rs_paths.cpp src/grid_a_star.cpp src/trailerlib.cpp src/trailer_hybrid_a_star.cpp)
target_link_libraries(lib_cpp_hybrid_a_star  Eigen3::Eigen ${PYTHON_LIBRARIES} Python3::NumPy)

add_executable(hybrid_a_star src/main.cpp)
target_link_libraries(hybrid_a_star  Eigen3::Eigen lib_cpp_hybrid_a_star  ${PYTHON_LIBRARIES} Python3::NumPy)

# find_package(IPOPT REQUIRED)
# include_directories(${IPOPT_INCLUDE_DIRS})

set(IPOPT_LIBRARY_DIRS /usr/local/lib)
set(IPOPT_INCLUDE_DIRS /usr/local/include/coin)


# Add Matplotlib C++ header
include_directories(${CMAKE_SOURCE_DIR} ${IPOPT_INCLUDE_DIRS})  # Assuming matplotlibcpp.h is in the source directory

# Add project source files
add_executable(TrajectoryOptimization src/try_ipopt.cpp)

# Link libraries
target_include_directories(TrajectoryOptimization PRIVATE ${CPPAD_INCLUDE_DIR} ${IPOPT_INCLUDE_DIRS})
target_link_libraries(TrajectoryOptimization ${IPOPT_LIBRARY_DIRS} Eigen3::Eigen Python3::NumPy ${PYTHON_LIBRARIES})

# Link Python for matplotlibcpp (assuming Python 3)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
target_include_directories(TrajectoryOptimization PRIVATE ${Python3_INCLUDE_DIRS})
target_link_libraries(TrajectoryOptimization ${Python3_LIBRARIES} ipopt )


