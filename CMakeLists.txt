cmake_minimum_required(VERSION 3.7.0...${CMAKE_VERSION})


if (BUILD_REMESH)
  add_compile_definitions(BUILD_REMESH)
endif()

#SET (BUILD_REMESH ON)
if (BUILD_OMEGA_H STREQUAL(ON))
message("SETTING OMEGA_H TO ON, FORCED REMESH TO ON")
else()
message("SETTING OMEGA_H TO OFF")
endif()

set (BUILD_MMG OFF)

if (BUILD_MMG)
  add_compile_definitions(MMG_BUILD)
endif()

if (BUILD_GPU)
  add_compile_definitions(CUDA_BUILD)
else()
  # include_directories(${PROJECT_SOURCE_DIR}/lib/Math/include)
endif (BUILD_GPU)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src/cuda)
include_directories(${PROJECT_SOURCE_DIR}/src/include)

MESSAGE ("CHECKING REMESH")


if (BUILD_MMG)
  MESSAGE ("ADDING MMG")
  add_subdirectory("lib/mmg-5.8.0")    
endif()


#AFTER REMESH
add_subdirectory (src)


include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src/include)

