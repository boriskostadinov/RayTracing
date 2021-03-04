cmake_minimum_required(VERSION 3.10)
project (bobi)

add_definitions(-D_CRT_SECURE_NO_WARNINGS)
add_definitions(-DNOMINMAX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(
	$ENV{CUDA_PATH}/include
)

link_directories(
	$ENV{CUDA_PATH}/lib/x64
)

set(HEADERS
)

include_directories(
	OpitX
	OpitX/include
	OpitX/sutil
	OpitX/sutil/tinyobjloader
	OpitX/sutil/rply-1.01
)
set(SOURCES
	"main.cpp"
	OpitX/sutil/Mesh.cpp
	OpitX/sutil/tinyobjloader/tiny_obj_loader.cc
)
set (C_SOURCES
	OpitX/sutil/rply-1.01/rply.c
)
#set_source_files_properties(OpitX/sutil/rply-1.01/rply.c PROPERTIES LANGUAGE CXX )

link_directories("${PROJECT_SOURCE_DIR}/OpitX/lib")

add_executable(${PROJECT_NAME} ${HEADERS} ${SOURCES} ${C_SOURCES})

target_link_libraries(${PROJECT_NAME} cuda.lib cudart.lib)
target_link_libraries(${PROJECT_NAME} optix.6.5.0.lib optix_prime.6.5.0.lib optixu.6.5.0.lib)