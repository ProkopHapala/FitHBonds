cmake_minimum_required(VERSION 3.10)

project(FitHBonds LANGUAGES C CXX)

option(RELEASE   "build optimized release (-O3/-Ofast) rather than debug? (-g -Og)"      OFF)
option(WITH_ASAN "use runtime memory sanitizer (asan) to trace segfaults and mem-leaks?" OFF)
option(WITH_OMP  "build with OpenMP paralelization for multiple CPUs"                  OFF)

message("OPTIONS: -DRELEASE=${RELEASE} -DWITH_ASAN=${WITH_ASAN} -DWITH_OMP=${WITH_OMP}")

# ==============================
# ==== FLAGS
# ==============================

if(RELEASE)
    set(CMAKE_BUILD_TYPE Release)
    add_definitions(-DRELEASEBUILD)
else()
    set(CMAKE_BUILD_TYPE Debug)
    add_definitions(-DDEBUGBUILD)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -fno-strict-aliasing")
endif()

message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

if(CMAKE_BUILD_TYPE MATCHES Debug)
    if(WITH_ASAN)
        set(AUX_COMPILE_FLAGS "-g -Og -march=native -fsanitize=address -fno-omit-frame-pointer")
    else()
        set(AUX_COMPILE_FLAGS "-g -Og -march=native")
    endif()
else()
    set(AUX_COMPILE_FLAGS "-w -Ofast -march=native -mtune=native")
endif()

if(WITH_OMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        set(AUX_COMPILE_FLAGS "${AUX_COMPILE_FLAGS} -fopenmp")
    else()
        message(WARNING "OpenMP not found, building without it")
    endif()
endif()

# Warning flags (copied from original; adjust as needed)
set(IGNORE_WARNING_FLAGS "-Wno-int-to-pointer-cast -Wno-switch -Wno-reorder -Wno-deprecated -Wno-maybe-uninitialized -Wno-char-subscripts -Wno-write-strings -Wno-format -Wno-parentheses -Wno-unused-but-set-variable -Wno-narrowing -Wno-unused-result -Wno-sign-compare -Wno-strict-aliasing -Wno-unused-variable -Wno-unused-value -Wno-comment -Wno-misleading-indentation")
set(WARNING_TO_ERROR "-Werror=return-type -Werror=init-self -Werror=uninitialized -Werror=implicit-fallthrough -Werror=tautological-compare -Werror=delete-non-virtual-dtor")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${AUX_COMPILE_FLAGS} ${IGNORE_WARNING_FLAGS} ${WARNING_TO_ERROR}")
message(STATUS "BUILD CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

# ==============================
# ==== DIRS
# ==============================

set(COMMON_RESOURCES "${CMAKE_SOURCE_DIR}/common_resources")

include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/math
    ${CMAKE_SOURCE_DIR}/molecular
)

# ==============================
# ==== BUILD TARGETS
# ==============================

# Shared library for FitREQ
add_library(FitREQ_lib SHARED libs/FitREQ_lib.cpp)
set_target_properties(FitREQ_lib PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/Build/libs/Molecular"
)

if(WITH_OMP AND OpenMP_CXX_FOUND)
    target_link_libraries(FitREQ_lib PUBLIC OpenMP::OpenMP_CXX)
endif()

# Optional: compile_commands.json for IDE support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    set(_CC_SRC "${CMAKE_BINARY_DIR}/compile_commands.json")
    set(_CC_DST "${CMAKE_SOURCE_DIR}/compile_commands.json")
    if(EXISTS ${_CC_DST} OR IS_SYMLINK ${_CC_DST})
        file(REMOVE ${_CC_DST})
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${_CC_SRC} ${_CC_DST})
endif()