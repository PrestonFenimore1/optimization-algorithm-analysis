cmake_minimum_required(VERSION 3.21)
project(source)

set(CMAKE_CXX_STANDARD 17)

add_executable(source main.cpp Population.cpp Population.h Equations.cpp Equations.h Algorithms.cpp Algorithms.h)
