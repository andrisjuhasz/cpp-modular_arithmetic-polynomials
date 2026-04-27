cmake_minimum_required(VERSION 3.20)

project(modular_arithmetic-polynomials LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(${PROJECT_NAME}
    src/demo.cpp
)

target_include_directories(${PROJECT_NAME}
    PRIVATE
        include
)

target_compile_options(${PROJECT_NAME}
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
)
