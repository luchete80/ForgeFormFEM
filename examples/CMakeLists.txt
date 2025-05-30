CMAKE_MINIMUM_REQUIRED (VERSION 3.18)

set (BUILD_EXAMPLES OFF)

project (ForgeFormFEM LANGUAGES CXX )

include_directories(${PROJECT_SOURCE_DIR})
include_directories(../src/include)
include_directories(${EIGEN_INC_DIR})
message ("EIGEN DIR" ${EIGEN_INC_DIR})



#LINUX
#find_package(PETSc REQUIRED)


add_executable(linear_test_matrix linear_test_matrix.cpp)
target_include_directories(linear_test_matrix PRIVATE ../src)

add_executable(test_1elem 
                    test_1elem.cpp
                    )
                    
target_link_libraries(test_1elem ForgeFormFEMlib)

#add_executable(linear_test_petsc linear_test_petsc.cpp)

#LINUX
#target_include_directories(linear_test_petsc PRIVATE ${PETSC_INCLUDES})
#target_link_libraries(linear_test_petsc PRIVATE ${PETSC_LIBRARIES})

add_executable(linear_test_eigen test_linear_eigen2_sparseLU.cpp)
