cmake_minimum_required(VERSION 3.14)

project(SimulationStaticLib LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(
    D:/Matlab/extern/include
    D:/Git/ThreeBodyFractal)


add_library(SimulationStaticLib STATIC
  SimulationStaticLib.cpp
  SimulationStaticLib.h
  ../defines.h
  ../Simulator.h
  ../Simulator.cpp
)

target_compile_definitions(SimulationStaticLib PRIVATE SIMULATIONSTATICLIB_LIBRARY)
