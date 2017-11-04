cmake_minimum_required(VERSION 2.6)
project(bobyqa)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Werror -std=c11 -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wold-style-cast -Werror -std=c++11 -fPIC")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)
include_directories(include)
set(SOURCES
    src/altmov.cpp
    src/bobyqa.cpp
    src/trsbox.cpp
    src/update.cpp
)
add_library(bobyqa_static STATIC ${SOURCES})
set_target_properties(bobyqa_static PROPERTIES OUTPUT_NAME bobyqa)
add_library(bobyqa_shared SHARED ${SOURCES})
set_target_properties(bobyqa_shared PROPERTIES OUTPUT_NAME bobyqa)
add_executable(example_function src/example_function.c)
target_link_libraries(example_function bobyqa_static)
add_executable(example_closure src/example_closure.cpp)
target_link_libraries(example_closure bobyqa_static)
enable_testing()
add_executable(test_bobyqa tests/bobyqa.cpp)
target_link_libraries(test_bobyqa bobyqa_static)
add_test(run_test_bobyqa bin/test_bobyqa)
install(TARGETS bobyqa_static ARCHIVE DESTINATION lib)
install(TARGETS bobyqa_shared LIBRARY DESTINATION lib)
install(FILES include/bobyqa.h DESTINATION include)
