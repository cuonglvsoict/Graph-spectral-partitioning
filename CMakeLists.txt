cmake_minimum_required(VERSION 3.17.5)

project(spectral_partition)

set(CMAKE_CXX_STANDARD 11)

file(GLOB ASSORTMENT_SRC
     "*.h"
     "*.cpp"
)

add_executable(${PROJECT_NAME} ${ASSORTMENT_SRC})