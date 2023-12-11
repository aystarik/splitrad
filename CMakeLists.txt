cmake_minimum_required(VERSION 2.8)
add_compile_options(-ffast-math)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project(splitrad)
add_executable(${PROJECT_NAME} "main.cpp")
add_executable("prune" "prune.cpp")
