cmake_minimum_required(VERSION 3.17.5)

project(assortment)

set(CMAKE_CXX_STANDARD 11)

file(GLOB ASSORTMENT_SRC
     "*.h"
     "*.cpp"
)

add_executable(${PROJECT_NAME} ${ASSORTMENT_SRC})

target_include_directories(${PROJECT_NAME} PRIVATE 
	${CMAKE_CURRENT_SOURCE_DIR}/../cplex-core/cplex/include
	${CMAKE_CURRENT_SOURCE_DIR}/../cplex-core/concert/include
)

target_link_libraries(${PROJECT_NAME}
	dl
	pthread
	${CMAKE_CURRENT_SOURCE_DIR}/../cplex-core/concert/lib/x86-64_linux/static_pic/libconcert.a 
	${CMAKE_CURRENT_SOURCE_DIR}/../cplex-core/cplex/lib/x86-64_linux/static_pic/libilocplex.a
	${CMAKE_CURRENT_SOURCE_DIR}/../cplex-core/cplex/lib/x86-64_linux/static_pic/libcplex.a
)