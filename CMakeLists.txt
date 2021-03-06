cmake_minimum_required (VERSION 2.6)

project (L1PulseAnalysis)
execute_process(COMMAND root-config --cflags OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE ROOTCFLAGS)
set(ROOT_LIBRARIES -L/usr/lib/i386-linux-gnu -lCore -lCint -lRIO -lNet -lHist -lGraf -lGpad -lTree -lRint -lMatrix -lMathCore -lThread -pthread -lm -ldl -rdynamic)
set(EIGEN_DIR /usr/include/eigen3)

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
add_executable (makeCaen1742Template ${PROJECT_SOURCE_DIR}/src/makeTemplateCaen1742.cxx)
add_executable (makeCaen5730Template ${PROJECT_SOURCE_DIR}/src/makeTemplateCaen5730.cxx)

target_link_libraries(pulseAnalysis projectlibs)
target_link_libraries(makeCaen1742Template projectlibs)
target_link_libraries(makeCaen5730Template projectlibs)

install(TARGETS makeCaen1742Template makeCaen5730Template pulseAnalysis DESTINATION ${PROJECT_SOURCE_DIR}/bin/)