cmake_minimum_required (VERSION 2.6)
project (bobyqa)
set (CMAKE_C_FLAGS "-Wall -lm -O3 -g")
set (CMAKE_CXX_FLAGS "-Wall -lm -O3 -g -std=c++11")
include_directories (include)
add_library (bobyqa_static STATIC src/bobyqa.cpp include/bobyqa.h)
set_target_properties (bobyqa_static PROPERTIES
    OUTPUT_NAME bobyqa
    ARCHIVE_OUTPUT_DIRECTORY lib)
add_library (bobyqa_shared SHARED src/bobyqa.cpp include/bobyqa.h)
set_target_properties (bobyqa_shared PROPERTIES
    OUTPUT_NAME bobyqa
    LIBRARY_OUTPUT_DIRECTORY lib)
add_executable (example_function src/example_function.c include/bobyqa.h)
set_target_properties (example_function PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY bin)
target_link_libraries (example_function bobyqa_static)
add_executable (example_closure src/example_closure.cpp include/bobyqa.h)
set_target_properties (example_closure PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY bin)
target_link_libraries (example_closure bobyqa_static)