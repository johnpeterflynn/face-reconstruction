cmake_minimum_required(VERSION 2.8)

project(face_recon_2019)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
add_definitions(-DPROJECT_DIR="${PROJECT_SOURCE_DIR}")

include_directories(include)
include_directories(thirdparty/Sophus)

# Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# Ceres Solver
find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})

find_package( OpenCV 3 QUIET COMPONENTS core features2d)
if (OpenCV_FOUND)
  message(STATUS "Found OpenCV 3 headers in: ${OpenCV_INCLUDE_DIRS}")
else()
  message(STATUS "Didn't find OpenCV 3, looking for 2")
  # HINTS needed for macos homebrew (install with `brew install opencv@2`)
  find_package( OpenCV 2 REQUIRED COMPONENTS core features2d
      HINTS /usr/local/opt/opencv@2 )
  message(STATUS "Found OpenCV 2 headers in: ${OpenCV_INCLUDE_DIRS}")
endif()

find_package(nabo)

set(SOURCES
    src/main.cc
    src/facemodel.cpp
    src/rgbdscan.cpp
    src/facesolver.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} OpenMeshCore ${OpenCV_LIBS} nabo ceres)
