CMAKE_MINIMUM_REQUIRED (VERSION 3.18)

set (BUILD_EXAMPLES OFF)

project (ForgeFormFEM LANGUAGES CXX )

include_directories(${PROJECT_SOURCE_DIR})
include_directories(../src/include)

#LINUX
find_package(PETSc REQUIRED)


add_executable(linear_test_matrix linear_test_matrix.cpp)

add_executable(linear_test_petsc linear_test_petsc.cpp)

#LINUX
target_include_directories(linear_test_petsc PRIVATE ${PETSC_INCLUDES})
target_link_libraries(linear_test_petsc PRIVATE ${PETSC_LIBRARIES})
