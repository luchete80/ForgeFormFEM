CMAKE_MINIMUM_REQUIRED (VERSION 3.18)

set (BUILD_EXAMPLES OFF)

project (ForgeFormFEM LANGUAGES CXX )

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/include)

include_directories(../lib/LSDynaReader/src/)

include_directories(../lib/omega_h-9.34.13/src) #FOR: Omega_h_config.h
set (MMG_ROOT "${CMAKE_BINARY_DIR}/lib/mmg-5.8.0/")
link_directories("$MMG_ROOT/lib")
include_directories("${CMAKE_BINARY_DIR}/lib/mmg-5.8.0/include") #FOR: Omega_h_config.h

set(OMEGA_H_DIR "${CMAKE_BINARY_DIR}/lib/omega_h-9.34.13/src")  # Adjust this path
link_directories(${OMEGA_H_DIR})

set(CMAKE_PREFIX_PATH "/path/to/omega_h")
include_directories(${CMAKE_BINARY_DIR}/lib/omega_h-9.34.13/src) #FOR: Omega_h_config.h


INCLUDE      (FindOpenMP)
message ("BUILD GPU set to OFF")
if (OPENMP_FOUND)
    message ("OPENMP FOUND!!!")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
message ("OPENMP NOT FOUND")
endif()
set_source_files_properties(*.C 
                          Matrices.C 
                           include/*.h src/include/defs.h src/include/utils.h
                          ${PROJECT_SOURCE_DIR}/src/include/*.cuh

                          PROPERTIES LANGUAGE CXX)
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3")



set (SOURCES                          
                         #Contact.C
                         Domain_d.C
                         #Matrices.C
                         #Mechanical.C
                         Mesh.C
                         Matrices.C
                         Solver.C
                         #input.C
                         ../lib/LSDynaReader/src/lsdynaReader.cpp
                         #VTKWriter.C
                   )

if (BUILD_REMESH)
 set (SOURCES "${SOURCES}" ReMesher.C
                           ReMesher_mmg.C)
endif()

                   
add_executable(ForgeFormFEM ForgeFormFEM.C
                         "${SOURCES}"
               )


