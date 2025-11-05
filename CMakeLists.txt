cmake_minimum_required(VERSION 3.16)
project(pileup-events CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(ENABLE_DEBUG_FLAGS "Enable extra debug compile options" OFF)

# direct paths to deps
set(CXXOPTS_ROOT "" CACHE PATH "Root of cxxopts (contains include/ or cxxopts.hpp)")
set(CXXOPTS_INCLUDE_DIR "" CACHE PATH "Path to directory that contains cxxopts.hpp")
set(HTSLIB_ROOT "" CACHE PATH "Root/prefix of htslib (contains include/ and lib/)")
set(HTSLIB_INCLUDE_DIR "" CACHE PATH "Path to directory that contains htslib/hts.h")
set(HTSLIB_LIBRARY "" CACHE FILEPATH "Full path to libhts.{so,dylib,a}")



add_executable(pileup-events
  src/main.cpp
)

# --- Dependencies: prefer pkg-config, otherwise search common macOS/Linux prefixes ---
find_package(PkgConfig QUIET)

# Common prefixes (Linux, macOS Intel/Apple Silicon, MacPorts)
set(_COMMON_HINT_PREFIXES
  /usr /usr/local /opt/local /opt/homebrew
)

# ---------------- cxxopts (header-only) ----------------
set(CXXOPTS_TARGET "")
set(_cxxopts_inc "")

# try user, then pkg-config, then general
if(CXXOPTS_INCLUDE_DIR)
  set(_cxxopts_inc "${CXXOPTS_INCLUDE_DIR}")
else()
  if(PkgConfig_FOUND)
    pkg_check_modules(CXXOPTS QUIET IMPORTED_TARGET cxxopts)
    if(CXXOPTS_FOUND)
      set(CXXOPTS_TARGET PkgConfig::CXXOPTS)
    endif()
  endif()
  find_path(_cxxopts_inc
    NAMES cxxopts.hpp
    HINTS ${CXXOPTS_ROOT} ${_COMMON_HINT_PREFIXES}
    PATH_SUFFIXES include
  )
endif()

if (NOT CXXOPTS_TARGET)
  add_library(cxxopts INTERFACE)
  add_library(cxxopts::cxxopts ALIAS cxxopts)
  target_include_directories(cxxopts INTERFACE "${_cxxopts_inc}")
  set(CXXOPTS_TARGET cxxopts::cxxopts)
endif()

if(NOT CXXOPTS_TARGET)
  message(FATAL_ERROR "Could not find cxxopts directly or via pkg-config. Set CXXOPTS_ROOT or CXXOPTS_INCLUDE_DIR.")
endif()


# ---------------- htslib (library + headers) ----------------
# --- htslib discovery (user -> pkg-config -> general) ---

set(HTSLIB_TARGET "")
set(_htslib_inc "")
set(_htslib_lib "")

if(HTSLIB_INCLUDE_DIR)
  set(_htslib_inc "${HTSLIB_INCLUDE_DIR}")
endif()
if(HTSLIB_LIBRARY)
  set(_htslib_lib "${HTSLIB_LIBRARY}")
endif()

if(NOT HTSLIB_TARGET AND PkgConfig_FOUND)
  # assuming the pkg-config module is called "htslib"
  pkg_check_modules(HTSLIB QUIET IMPORTED_TARGET htslib)
  if(HTSLIB_FOUND)
    set(HTSLIB_TARGET PkgConfig::HTSLIB)
  endif()
endif()

if(NOT HTSLIB_TARGET)
  if(NOT _htslib_inc)
    find_path(_htslib_inc
      NAMES htslib/hts.h
      HINTS ${HTSLIB_ROOT} ${_COMMON_HINT_PREFIXES}
      PATH_SUFFIXES include
    )
  endif()

  if(NOT _htslib_lib)
    find_library(_htslib_lib
      NAMES htslib
      HINTS ${HTSLIB_ROOT} ${_COMMON_HINT_PREFIXES}
      PATH_SUFFIXES lib lib64
    )
  endif()
endif()

if(NOT HTSLIB_TARGET AND _htslib_inc AND _htslib_lib)
  add_library(htslib IMPORTED GLOBAL)
  add_library(htslib::htslib ALIAS htslib)
  set_target_properties(htslib PROPERTIES
    IMPORTED_LOCATION "${_htslib_lib}"
    INTERFACE_INCLUDE_DIRECTORIES "${_htslib_inc}"
  )
  set(HTSLIB_TARGET htslib::htslib)
endif()

if(NOT HTSLIB_TARGET)
  message(FATAL_ERROR
    "Could not find htslib directly or via pkg-config. "
    "Set HTSLIB_ROOT or both HTSLIB_INCLUDE_DIR and HTSLIB_LIBRARY."
  )
endif()
# --- end deps ---

target_compile_options(pileup-events PRIVATE
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

if(ENABLE_DEBUG_FLAGS)
  target_compile_options(pileup-events PRIVATE
    -fsanitize=address,undefined
    -fno-omit-frame-pointer
    -O0 -g
  )
  target_link_options(pileup-events PRIVATE
    -fsanitize=address,undefined
  )
endif()

target_link_libraries(pileup-events PRIVATE
  ${CXXOPTS_TARGET}
  ${HTSLIB_TARGET}
)

