cmake_minimum_required(VERSION 2.8)
project(pLLL)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/)

#include(CheckLanguage)
#check_language(CXX)
#if(NOT CMAKE_CXX_COMPILER)
#    message(FATAL_ERROR "No C++ compiler found!")
#endif(NOT CMAKE_CXX_COMPILER)

if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Use -DCMAKE_BUILD_TYPE=Release or -DCMAKE_BUILD_TYPE=Debug to set a specific build type. Now using default compile settings.")
endif(NOT CMAKE_BUILD_TYPE)
if(CMAKE_BUILD_TYPE STREQUAL Debug)
    add_definitions(-DDEBUG)
endif(CMAKE_BUILD_TYPE STREQUAL Debug)

if(FASTCOMPILE)
    set(PLLL_CONFIG_NO_ARITHMETIC_RATIONAL True)
    set(PLLL_CONFIG_NO_ARITHMETIC_REAL True)
    set(PLLL_CONFIG_NO_ARITHMETIC_DOUBLE True)
    #set(PLLL_CONFIG_NO_ARITHMETIC_LONGDOUBLE False)
    set(PLLL_CONFIG_NO_ARITHMETIC_DOUBLEDOUBLE True)
    set(PLLL_CONFIG_NO_ARITHMETIC_QUADDOUBLE True)
    
    #set(PLLL_CONFIG_NO_ARITHMETIC_BIGINT False)
    #set(PLLL_CONFIG_NO_ARITHMETIC_LONG False)
    
    set(PLLL_CONFIG_NO_GS_CLASSIC True)
    set(PLLL_CONFIG_NO_GS_CLASSICINT True)
    set(PLLL_CONFIG_NO_GS_GIVENS True)
    #set(PLLL_CONFIG_NO_GS_NUMSTABLE False)
    message(STATUS "Doing fast compile")
else(FASTCOMPILE)
    set(FASTCOMPILE False)
endif(FASTCOMPILE)

if(CPP11)
    message(STATUS "Trying to enable C++11 compilation.")
    set(PLLL_CPP11_COMPILATION True)
endif(CPP11)

include(CompilerOptions.txt)
include(Libraries.txt)

configure_file(${PROJECT_SOURCE_DIR}/plll/include/plll/config.hpp.in ${PROJECT_BINARY_DIR}/plll/include/plll/config.hpp @ONLY)

add_subdirectory(plll)
add_subdirectory(tests)
add_subdirectory(tools)
add_subdirectory(doc)
add_subdirectory(examples)

# uninstall target; see http://www.itk.org/Wiki/CMake_FAQ#Can_I_do_.22make_uninstall.22_with_CMake.3F
#configure_file(
#    "${PROJECT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
#    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
#    IMMEDIATE @ONLY)
#add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

message(STATUS "-----------------------------------------------------------------------------------------")
message(STATUS "Source DIR is ${PROJECT_SOURCE_DIR}")
message(STATUS "Binary DIR is ${PROJECT_BINARY_DIR}")
message(STATUS "Build type is ${CMAKE_BUILD_TYPE}")
message(STATUS "Fast compile mode is ${FASTCOMPILE}")
message(STATUS "Libs extension is '${SUFFIXLIB}'. Bins extension is '${SUFFIXBIN}'." )
message(STATUS "Installation prefix directory is " ${CMAKE_INSTALL_PREFIX})
message(STATUS "Host system is " ${CMAKE_HOST_SYSTEM} " with processor " ${CMAKE_HOST_SYSTEM_PROCESSOR})
message(STATUS "Target system is " ${CMAKE_SYSTEM} " with processor " ${CMAKE_SYSTEM_PROCESSOR})
