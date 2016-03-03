cmake_minimum_required (VERSION 2.6)

project (L1PulseAnalysis)
find_package(ROOT REQUIRED)
set(EIGEN_DIR /usr/local/include/eigen3)

add_definitions(-std=c++0x -Wall -Wextra ${ROOTCFLAGS})

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
include_directories(${PROJECT_SOURCE_DIR}/../italian-testbeam-daq/fast/core/include)

add_library(projectlibs STATIC ${libsrcs} ${utility})
target_link_libraries(projectlibs ${ROOT_LIBRARIES})

add_executable (pulseAnalysis ${PROJECT_SOURCE_DIR}/src/pulseAnalyzer.cxx)
add_executable (makeCaen5742Template ${PROJECT_SOURCE_DIR}/src/makeTemplateCaen5742.cxx)
add_executable (makeCaen5730Template ${PROJECT_SOURCE_DIR}/src/makeTemplateCaen5730.cxx)

target_link_libraries(pulseAnalysis projectlibs)
target_link_libraries(makeCaen5742Template projectlibs)
target_link_libraries(makeCaen5730Template projectlibs)

install(TARGETS makeCaen5742Template makeCaen5730Template pulseAnalysis DESTINATION ${PROJECT_SOURCE_DIR}/../bin/)