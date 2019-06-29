cmake_minimum_required(VERSION 2.8)

project(face_recon_2019)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
add_definitions(-DPROJECT_DIR="${PROJECT_SOURCE_DIR}")

# Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# Ceres Solver
#find_package(Ceres REQUIRED)
#include_directories(${CERES_INCLUDE_DIRS})

add_executable(${PROJECT_NAME} "main.cc")
target_link_libraries(${PROJECT_NAME} OpenMeshCore)
