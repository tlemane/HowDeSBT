project(HowDeSBT)

cmake_minimum_required(VERSION 3.10.0)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

if (NOT APPLE)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

find_package(ZLIB REQUIRED)

include(ExternalProject)
ExternalProject_add(CROAR
    PREFIX CROAR
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/thirdparty/CRoaring
    CMAKE_ARGS -DROARING_BUILD_STATIC=ON
    INSTALL_COMMAND "")

link_directories(${PROJECT_BINARY_DIR}/CROAR/src/CROAR-build/src)
include_directories(${PROJECT_SOURCE_DIR}/thirdparty/CRoaring/include/ ${ZLIB_INCLUDE_DIRS})

if (NOT SDSLLIB OR NOT SDSLINC OR NOT KMTRICKS_INC OR NOT KMTRICKS_LIB OR NOT LZ4INC OR NOT LZ4LIB)
    message(FATAL_ERROR "-DSDSLLIB=path -DSDSLINC=path -DKMTRICKS_INC=path -DKMTRICKS_LIB=path are required")
else()

    if(NOT OUTPUT)
        set(OUTPUT ${PROJECT_SOURCE_DIR}/bin)
    endif()

    include_directories(${SDSLINC} ${KMTRICKS_INC} ${LZ4INC})
    link_directories(${SDSLLIB} ${KMTRICKS_LIB} ${LZ4LIB})

    install(DIRECTORY DESTINATION ${OUTPUT})
    set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OUTPUT})

    file(GLOB_RECURSE ProjectFiles *.cc)

    set(LIBS kmtricks roaring sdsl pthread ZLIB::ZLIB lz4)
    set(DEPS CROAR SDSL)
    include_directories(${PROJECT_SOURCE_DIR})

    add_executable(km_howdesbt ${ProjectFiles})
    add_dependencies(km_howdesbt CROAR)
    target_link_libraries(km_howdesbt ${LIBS})

endif()




