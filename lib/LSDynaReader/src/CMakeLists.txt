PROJECT(test)
cmake_minimum_required(VERSION 2.6)
enable_language(CXX)
# Setting the compilers
# set (CMAKE_Fortran_COMPILER /usr/bin/gfortran)
# set (CMAKE_CXX_COMPILER /usr/bin/g++)
# Setting the flags
# set (CMAKE_CXX_FLAGS "-lgfortran")
# set_source_files_properties(main.f90 func.cpp PROPERTIES COMPILE_FLAGS -c)
# Making the executable
ADD_EXECUTABLE(reader         main_c.cpp lsdynaReader.cpp)
ADD_LIBRARY   (lsdynareader   main_c.cpp lsdynaReader.cpp)
# TARGET_LINK_LIBRARIES (test gfortran quadmath)
# SET_TARGET_PROPERTIES(test PROPERTIES LINKER_LANGUAGE Fortran)

# ADD_EXECUTABLE(reader main_reader.cpp NastranReader.cpp )

