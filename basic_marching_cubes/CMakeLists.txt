cmake_minimum_required(VERSION 3.8)
project(basic_marching_cubes)

set(CMAKE_CXX_STANDARD 17)

add_executable(basic_marching_cubes
			${PROJECT_SOURCE_DIR}/main.cpp ${PROJECT_SOURCE_DIR}/marching_cubes.cpp)