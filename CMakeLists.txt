cmake_minimum_required(VERSION 2.8.9)
project(cheater_recognition)


SET(CMAKE_INSTALL_PREFIX ./)

file(GLOB SOURCES "./src/*.cpp")

find_library(GSL_LIB gsl)
find_library(GSLCBLAS_LIB gslcblas)

add_executable(cheater_recognition ${SOURCES})

target_link_libraries(cheater_recognition PRIVATE "${GSL_LIB}")
target_link_libraries(cheater_recognition PRIVATE "${GSLCBLAS_LIB}")

install(TARGETS cheater_recognition DESTINATION bin)
