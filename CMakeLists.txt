cmake_minimum_required(VERSION 3.16)
project(spars_gmres_ilu)
enable_language(Fortran)

set(GCC_COVERAGE_COMPILE_FLAGS "-O3 -march=native")
set(GCC_COVERAGE_COMPILE_FLAGS "-g -fsanitize=address")
set(CMAKE_Fortran_FLAGS "${GCC_COVERAGE_COMPILE_FLAGS} -lblas")

find_package(BLAS REQUIRED)

add_executable(spars_gmres_ilu main.f90)
target_link_libraries(spars_gmres_ilu ${CMAKE_SOURCE_DIR}/libskit.a)
