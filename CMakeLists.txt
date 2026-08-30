cmake_minimum_required(VERSION 3.16)
project(ControlToolBox CXX)

if(MSVC)
    add_compile_options(/W4 /permissive-)
else()
    add_compile_options(-Wall -Wextra -Wshadow -Wreorder -Wunreachable-code)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_library(control_toolbox)

target_sources(control_toolbox PRIVATE
    src/matrix.cpp
    src/active_set.cpp 
    src/simplex.cpp
    src/mpc.cpp
    src/pid.cpp
    src/kalman_filter.cpp
)

target_include_directories(control_toolbox PUBLIC include)


add_executable(simplex_test tests/test_simplex.cpp)
target_link_libraries(simplex_test PRIVATE control_toolbox)

add_executable(mpc_test tests/mpc_test.cpp)
target_link_libraries(mpc_test PRIVATE control_toolbox)

add_executable(active_set_test tests/test_active_set.cpp)
target_link_libraries(active_set_test PRIVATE control_toolbox)

add_executable(pid_test tests/pid_test.cpp)
target_link_libraries(pid_test PRIVATE control_toolbox)

add_executable(kalman_demo tests/kalman_demo.cpp)
target_link_libraries(kalman_demo PRIVATE control_toolbox)