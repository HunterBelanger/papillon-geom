cmake_minimum_required(VERSION 3.11)

project(PapillonGeom
  VERSION "0.1.0"
  DESCRIPTION "Library for Constructive Solid Geometry modeling."
  LANGUAGES CXX
)

include(CMakePackageConfigHelpers)
include(FetchContent)

# Add options
option(PGEOM_SHARED "Build PapillonGeom as a shared library" ON)
option(PGEOM_INSTALL "Install the PapillonGeom library and header files" ON)
option(PGEOM_TESTS "Build PapillonGeom tests" OFF)

# List of source files for PapillonGeom
set(PGEOM_SOURCE_LIST src/version.cpp)

# Only let a static library be built on Windows
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  if(PGEOM_SHARED)
    message("Only a static library may be built for Windows.")
    message("Setting PGEOM_SHARED to OFF.")
    set(PGEOM_SHARED OFF)
  endif()
endif()

# Create PapillonGeom library target
if(PGEOM_SHARED)
  add_library(PapillonGeom SHARED ${PGEOM_SOURCE_LIST})
else()
  add_library(PapillonGeom STATIC ${PGEOM_SOURCE_LIST})
endif()

# Add alias to make more friendly with FetchConent
add_library(PapillonGeom::PapillonGeom ALIAS PapillonGeom)

# Add public interface includes directory
target_include_directories(PapillonGeom PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# Require C++20 standard
target_compile_features(PapillonGeom PUBLIC cxx_std_20)


if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC") # Comile options for Windows
  target_compile_options(PapillonGeom PRIVATE /W4)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU") # Compile options for GCC
  target_compile_options(PapillonGeom PRIVATE -W -Wall -Wextra -Wpedantic)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang") # Compile options for Clang
  target_compile_options(PapillonGeom PRIVATE -W -Wall -Wextra -Wpedantic)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel") # Compile options for Intel
  target_compile_options(PapillonGeom PRIVATE -W -Wall -Wextra -Wpedantic)
endif()

# If building tests, add the tests subdirectory
if(PGEOM_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

# Install PapillonGeom
if(PGEOM_INSTALL)
  configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/PapillonGeomConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/PapillonGeomConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PapillonGeom
  )

  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/PapillonGeomConfigVersion.cmake"
    VERSION "${PapillonGeom_VERSION}"
    COMPATIBILITY SameMinorVersion
  )

  install(TARGETS PapillonGeom
          EXPORT PapillonGeomTargets
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
          INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/PapillonGeom DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

  install(EXPORT PapillonGeomTargets
          FILE PapillonGeomTargets.cmake
          NAMESPACE PapillonGeom::
          DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PapillonGeom
  )

  install(FILES
            "${CMAKE_CURRENT_BINARY_DIR}/PapillonGeomConfig.cmake"
            "${CMAKE_CURRENT_BINARY_DIR}/PapillonGeomConfigVersion.cmake"
          DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PapillonGeom
  )
endif()
