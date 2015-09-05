cmake_minimum_required (VERSION 2.6)
project (bobyqa)
set (CMAKE_C_FLAGS "-Wall -lm -O3 -g")
set (CMAKE_CXX_FLAGS "-Wall -lm -O3 -g")
include_directories (include)
add_library (bobyqa_static STATIC src/bobyqa.cpp include/bobyqa.h)
set_target_properties (bobyqa_static PROPERTIES
    OUTPUT_NAME bobyqa
    ARCHIVE_OUTPUT_DIRECTORY lib)
add_library (bobyqa_shared SHARED src/bobyqa.cpp include/bobyqa.h)
set_target_properties (bobyqa_shared PROPERTIES
    OUTPUT_NAME bobyqa
    LIBRARY_OUTPUT_DIRECTORY lib)
add_executable (example src/example.c include/bobyqa.h)
set_target_properties (example PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY bin)
target_link_libraries (example bobyqa_static)
