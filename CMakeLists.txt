cmake_minimum_required(VERSION 3.10)

project(EDMFitter)

#include debug hooks:
#SET(CMAKE_BUILD_TYPE "Debug")

SET(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/install)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

#message(STATUS "Source directory ${PROJECT_SOURCE_DIR}")
#message(STATUS "Binary directory ${PROJECT_BINARY_DIR}")
#message(STATUS "Install directory ${CMAKE_INSTALL_PREFIX}")
#message(STATUS "Header directory ${CMAKE_INSTALL_PREFIX}/include")
#message(STATUS "Library directory ${CMAKE_INSTALL_PREFIX}/lib")

##ROOT stuff:
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS)
#include_directories(${ROOT_INCLUDE_DIRS})
#add_definitions(${ROOT_CXX_FLAGS})

#my header files:
#include_directories(include)

#all the source code:
add_subdirectory(3rdparty)
add_subdirectory(src)


