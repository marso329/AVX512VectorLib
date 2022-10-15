cmake_minimum_required (VERSION 3.9.1)

project (AVX512VectorLib VERSION 0.0.1 DESCRIPTION "AVX512VectorLib project")

set(CMAKE_BUILD_TYPE Release)

include(GNUInstallDirs)


#All gcc/icc flags, feel free to add flags if necessary
set (CMAKE_CXX_FLAGS "-std=c++17 -g -march=skylake-avx512 ")

#All warning flags, yes we are pedantic
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra   -Wcast-align   -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-include-dirs -Wnoexcept   -Wredundant-decls  -Wsign-promo -Wstrict-null-sentinel   -Wno-unused -Wno-variadic-macros -Wno-parentheses -fdiagnostics-show-option -Wno-deprecated-declarations")

#All libraries that must be linked
set (CMAKE_SHARED_LINKER_FLAGS "-lrt -lpthread -ldl -lstdc++fs")

#I like the output files to end up in .lib
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)


add_library(AVX512VectorLib STATIC
    src/AVX512VectorLib.cxx
)

set_target_properties(AVX512VectorLib PROPERTIES VERSION ${PROJECT_VERSION})

#Not needed but libTemplate.so.1 will symlink to libTemplate.so.0.0.1
set_target_properties(AVX512VectorLib PROPERTIES SOVERSION 1)

set_target_properties(AVX512VectorLib PROPERTIES PUBLIC_HEADER "include/AVX512VectorLib.hpp")

#include directory
target_include_directories(AVX512VectorLib PUBLIC include )

#Src directory
target_include_directories(AVX512VectorLib PUBLIC src)


#Rules on how to install it
install(TARGETS AVX512VectorLib
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_subdirectory(tests)
