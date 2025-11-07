cmake_minimum_required(VERSION 3.16)
project(pileup-events LANGUAGES CXX VERSION 0.0.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(ENABLE_DEBUG_FLAGS "Enable extra debug compile options" OFF)
option(MAKE_PY_BINDS "Build python bindings via SWIG" OFF)

# user paths to htslib if needed
set(HTSLIB_INCLUDE_DIR "" CACHE PATH "Path to directory that contains htslib/hts.h")
set(HTSLIB_LIBRARY "" CACHE FILEPATH "Absolute path to libhts.{so,dylib,a}")

# //----- find deps ---//
find_package(PkgConfig QUIET)

# ------- cxxopts -------
# download directly since light, header only
include(FetchContent)

message(STATUS "Fetching cxxopts dependency from remote.")
FetchContent_Declare(
  cxxopts
  GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
  GIT_TAG        v3.3.1  # what I developed against, but any 3.x.x should be fine
)
FetchContent_MakeAvailable(cxxopts)

# ------- htslib -------
set(HTSLIB_TARGET "")
set(_htslib_inc "")
set(_htslib_lib "")

if(HTSLIB_INCLUDE_DIR AND HTSLIB_LIBRARY)
  message(STATUS "Using user-specified htslib install. Must be >=1.14, but this is not enforced.")
  set(_htslib_inc "${HTSLIB_INCLUDE_DIR}")
  set(_htslib_lib "${HTSLIB_LIBRARY}")
endif()

if(_htslib_inc AND _htslib_lib)
  add_library(htslib UNKNOWN IMPORTED GLOBAL)
  add_library(htslib::htslib ALIAS htslib)
  set_target_properties(htslib PROPERTIES
    IMPORTED_LOCATION "${_htslib_lib}"
    INTERFACE_INCLUDE_DIRECTORIES "${_htslib_inc}"
  )
  set(HTSLIB_TARGET htslib::htslib)
endif()

if(NOT HTSLIB_TARGET AND PkgConfig_FOUND)
  # assuming the pkg-config module is called "htslib"
  pkg_check_modules(HTSLIB QUIET IMPORTED_TARGET htslib>=1.14)
  if(HTSLIB_FOUND)
    message(STATUS "Located htslib>=1.14 install via pkg-config.")
    set(HTSLIB_TARGET PkgConfig::HTSLIB)
  endif()
endif()

if(NOT HTSLIB_TARGET)
  message(FATAL_ERROR
    "Could not find htslib directly or via pkg-config. "
    "Set both HTSLIB_INCLUDE_DIR and HTSLIB_LIBRARY to use user install."
  )
endif()
# //--- end deps ---//

add_library(pev_core INTERFACE)
target_include_directories(pev_core INTERFACE include)

if(ENABLE_DEBUG_FLAGS)
  target_compile_options(pev_core INTERFACE
      -Wall
      -Wextra
      -Wpedantic
      -Wshadow
      -Wconversion
      -Wsign-conversion
      -Wfloat-equal
      -Wnon-virtual-dtor
      -Wold-style-cast
      -Woverloaded-virtual
      -Wnull-dereference
      -Wdouble-promotion
      -Wformat=2
  )
  target_compile_options(pev_core PRIVATE
    -fsanitize=address,undefined
    -fno-omit-frame-pointer
    -O0 -g
  )
  target_link_options(pev_core PRIVATE
    -fsanitize=address,undefined
  )
endif()

target_link_libraries(pev_core INTERFACE
  ${HTSLIB_TARGET}
)

set_target_properties(pev_core PROPERTIES POSITION_INDEPENDENT_CODE ON)


add_executable(pileup-events
  src/main.cpp
)
target_link_libraries(pileup-events
  PRIVATE
  cxxopts::cxxopts
  pev_core
)


# SWIG
if(MAKE_PY_BINDS)
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
  find_package(SWIG 4.0 REQUIRED)
  include(UseSWIG)
  set(CMAKE_SWIG_FLAGS -c++)

  set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/swig/pev_bindings.i
    PROPERTIES CPLUSPLUS ON)

  set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR}/python)
  swig_add_library(pev_bindings
    TYPE MODULE
    LANGUAGE python
    SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/swig/pev_bindings.i
  )

  set_target_properties(pev_bindings PROPERTIES
    SWIG_USE_TARGET_INCLUDE_DIRECTORIES ON
    OUTPUT_NAME "_pev"                        # import pev_bindings; loads _pev
  )
  target_include_directories(pev_bindings PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src              # if wrap.hpp is beside sources
  )
  target_link_libraries(pev_bindings PRIVATE pev_core Python3::Module)
endif(MAKE_PY_BINDS)


