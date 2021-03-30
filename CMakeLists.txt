cmake_minimum_required(VERSION 3.5)

project("Iximiel's Minimization Hub" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/include)
add_subdirectory(src)

