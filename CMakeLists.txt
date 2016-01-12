cmake_minimum_required (VERSION 2.6)

project (L1PulseAnalysis)
set(CMAKE_MODULE_PATH /home/newg2/Packages/root/etc/cmake/)
find_package(ROOT REQUIRED)
set(EIGEN_DIR /usr/include/eigen3)

add_definitions(-std=c++0x -Wall)

file(GLOB libsrcs 
	${PROJECT_SOURCE_DIR}/templateFitter/src/*.cxx
	${PROJECT_SOURCE_DIR}/json11/json11.cpp)
set(utility ${PROJECT_SOURCE_DIR}/src/utility.cxx)
set(projectincludes  ${PROJECT_SOURCE_DIR}/include 
	${PROJECT_SOURCE_DIR}/templateFitter/src/ 
	${PROJECT_SOURCE_DIR}/json11)

include_directories(${projectincludes})
include_directories(${ROOT_INCLUDE_DIRS}) 
include_directories(${EIGEN_DIR})

add_library(projectlibs STATIC ${libsrcs} ${utility})
target_link_libraries(projectlibs ${ROOT_LIBRARIES})

add_executable (pulseAnalysis ${PROJECT_SOURCE_DIR}/src/pulseAnalyzer.cxx)
add_executable (makeCaenTemplate ${PROJECT_SOURCE_DIR}/src/makeTemplateCaen1742.cxx)

target_link_libraries(pulseAnalysis projectlibs)
target_link_libraries(makeCaenTemplate projectlibs)

install(TARGETS makeCaenTemplate pulseAnalysis DESTINATION ${PROJECT_SOURCE_DIR}/../bin/)