cmake_minimum_required (VERSION 3.5)
project (lib_cpp_hybrid_a_star)

set(CMAKE_BUILD_TYPE "DEBUG")
set(CMAKE_CXX_FLAGS "-std=c++17")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g --debug-cpp --trace-expand")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_MATH_DISABLE_FLOAT128 -O1 -Wall -g --debug")


find_package(Eigen3 REQUIRED)

include_directories( 
    SYSTEM 
    include
    ${Eigen3_INCLUDE_DIRS}
)

add_library(lib_cpp_hybrid_a_star src/rs_paths.cpp)
target_link_libraries(lib_cpp_hybrid_a_star  Eigen3::Eigen)

add_executable(hybrid_a_star src/main.cpp)
target_link_libraries(hybrid_a_star  Eigen3::Eigen lib_cpp_hybrid_a_star)
