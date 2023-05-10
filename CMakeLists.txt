cmake_minimum_required(VERSION 3.18)

project(PyTrackBead VERSION 0.1)
set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_EXTENSIONS ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_DEBUG "-g")
add_compile_options(-ggdb3)

# Common installation directories
include(GNUInstallDirs)

# Use -fPIC even if statically compiled
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake
                      ${CMAKE_MODULE_PATH})


find_package(Python3 COMPONENTS Interpreter Development)
find_package(pybind11 REQUIRED)
find_package(Python3 COMPONENTS Interpreter)
find_package(Threads REQUIRED)






###################################################
# zwo_asi: library for higher level calls to the #
# ZWO library                                     #
###################################################


  add_library(Tracker
    src/tracking/Bead.cpp
    src/tracking/fftbtl32n.c
    src/tracking/fillibbt.c
    src/tracking/Tracker.cpp
    src/tracking/track_util.c
  )



 target_include_directories(Tracker PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/tracking>
   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
   $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
 )

 target_link_libraries(Tracker Threads::Threads m )

#################################
# python bindings over zwo_asi #
#################################

if(CALL_FROM_SETUP_PY)
    set(BEAD_TRACKER_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
else()
    set(BEAD_TRACKER_INSTALL_PREFIX ${Python3_SITELIB})
endif()
pybind11_add_module(pybind11_bindings MODULE srcpy/bindings.cpp)
set_target_properties(pybind11_bindings PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bead_tracker
  OUTPUT_NAME "bindings"
  )
target_link_libraries(pybind11_bindings
  PUBLIC Tracker
  )

###########
# install #
###########

install(
    TARGETS  Tracker #MODIFIED MR
    EXPORT TrackerTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

install(
    TARGETS pybind11_bindings
    COMPONENT bindings
    LIBRARY DESTINATION ${BEAD_TRACKER_INSTALL_PREFIX}
    ARCHIVE DESTINATION ${BEAD_TRACKER_INSTALL_PREFIX}
    RUNTIME DESTINATION ${BEAD_TRACKER_INSTALL_PREFIX}
    )
