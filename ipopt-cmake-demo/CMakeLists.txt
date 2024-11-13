
cmake_minimum_required(VERSION 2.8)
project(ipopt_example)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g")
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

set(IPOPT_LIBRARY_DIRS /usr/local/lib)
set(IPOPT_INCLUDE_DIRS /usr/local/include/coin)


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


message("ipopt binary path is: " ${IPOPT_LIBRARY_DIRS} )
message("ipopt include path is: " ${IPOPT_INCLUDE_DIRS} )

include_directories(${IPOPT_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)
link_directories(${IPOPT_LIBRARY_DIRS})

add_library(ipopt_interface1 ipopt_interface.cpp)
add_library(ipopt_interface ipopt_interface_traj.cpp)
target_link_libraries(ipopt_interface Eigen3::Eigen ${PYTHON_LIBRARIES} Python3::NumPy)
target_link_libraries(ipopt_interface1 Eigen3::Eigen ${PYTHON_LIBRARIES} Python3::NumPy)

add_executable(ipopt_example main.cpp)
add_executable(ipopt_orign main_original.cpp)

target_link_libraries(ipopt_orign  ipopt_interface1)
target_link_libraries(ipopt_example ipopt_interface ipopt Eigen3::Eigen ${PYTHON_LIBRARIES} Python3::NumPy)
target_link_libraries(ipopt_orign  ipopt Eigen3::Eigen ${PYTHON_LIBRARIES} Python3::NumPy)

