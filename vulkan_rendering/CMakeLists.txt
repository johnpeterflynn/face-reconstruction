cmake_minimum_required(VERSION 3.10)
project(vulkan_rendering)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11  " )
set(WORKING_DIRECTORY "E:\\Lex\\GameDev\\Projects\\3D Scanning Project\\vulkan_rendering")

#set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
#set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
#set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

add_executable(${PROJECT_NAME} main.cpp)

include_directories(libs)
include_directories(models)
include_directories(shaders)
include_directories(textures)

#include_directories(C:\\Users\\Lex\\Documents\\Visual Studio 2017\\Libraries\\glfw-3.3)

find_package(Vulkan REQUIRED)
#find_package(glfw3 3.3 REQUIRED)

if (VULKAN_FOUND)
    message(STATUS "Found Vulkan, Including and Linking now")
    include_directories(${Vulkan_INCLUDE_DIRS})
    target_link_libraries (${PROJECT_NAME} ${Vulkan_LIBRARIES} glfw)
endif (VULKAN_FOUND)
