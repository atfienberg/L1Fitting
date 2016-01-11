cmake_minimum_required (VERSION 2.6)

project (L1PulseAnalysis)
set(CMAKE_MODULE_PATH /home/newg2/Packages/root/etc/cmake/)
find_package(ROOT REQUIRED)
set(EIGEN_DIR /usr/include/eigen3)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cxx 
	${PROJECT_SOURCE_DIR}/templateFitter/src/*.cxx
	${PROJECT_SOURCE_DIR}/json11/json11.cpp)
set(projectincludes  ${PROJECT_SOURCE_DIR}/include 
	${PROJECT_SOURCE_DIR}/templateFitter/src/ 
	${PROJECT_SOURCE_DIR}/json11)

include_directories(${projectincludes})
include_directories(${ROOT_INCLUDE_DIRS}) 
include_directories(${EIGEN_DIR})

add_executable (pulseAnalysis ${sources})

target_link_libraries(pulseAnalysis ${ROOT_LIBRARIES})