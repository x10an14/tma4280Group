project(ex6)

cmake_minimum_required(VERSION 2.6)

# We use the C language
enable_language(C)
enable_language(Fortran)

# CMake Options to enable/disable MPI and OpenMP
option(ENABLE_OPENMP "Enable OpenMP support?" ON)
option(ENABLE_MPI    "Enable MPI support?" ON)

if(ENABLE_MPI)
  find_package(MPI)
endif()

if(MPI_FOUND)
  add_definitions(-DHAVE_MPI=1)
  set(INCLUDES ${INCLUDES} ${MPI_INCLUDE_PATH})
  set(DEPLIBS_C ${DEPLIBS_C} ${MPI_C_LIBRARIES})
  set(DEPLIBS_F ${DEPLIBS_F} ${MPI_Fortran_LIBRARIES})
endif()

if(ENABLE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    add_definitions(-DHAVE_OPENMP=1)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  endif()
endif()

add_definitions(${COMMON_DEFINITIONS})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMMON_C_FLAGS} -std=c99 -Wall")
include_directories(${COMMON_INCLUDE_PATHS})

# Send the include directories to the compiler
include_directories(${INCLUDES})

# Add a library with common sources //Not sure if we need common. Should try removing it
add_library(common ../fst.f)

# The C based ex6 application
add_executable(ex6 ex6.c)

## The Fortran based poisson application
#add_executable(poisson-f ../poisson.f)

# Link to libraries
target_link_libraries(ex6 common ${DEPLIBS_C})

# We want to link this as a C program
set_target_properties(ex6 PROPERTIES LINKER_LANGUAGE C)

# Link to libraries
target_link_libraries(common ${DEPLIBS_F})
