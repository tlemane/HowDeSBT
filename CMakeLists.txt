project(HowDeSBT)

cmake_minimum_required(VERSION 3.10.0)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

include(ExternalProject)
ExternalProject_add(CROAR
    PREFIX CROAR
    GIT_REPOSITORY https://github.com/RoaringBitmap/CRoaring
    INSTALL_COMMAND "")

link_directories(${PROJECT_BINARY_DIR}/CROAR/src/CROAR-build/)
include_directories(${PROJECT_BINARY_DIR}/CROAR/src/CROAR/include/)

if (NOT SDSLLIB OR NOT SDSLINC)
    ExternalProject_add(SDSL
        PREFIX SDSL
        GIT_REPOSITORY https://github.com/simongog/sdsl-lite
        INSTALL_COMMAND "")
    include_directories(${PROJECT_BINARY_DIR}/SDSL/src/SDSL-build/include)
    link_directories(${PROJECT_BINARY_DIR}/SDSL/src/SDSL-build/lib)
else()
    include_directories(SDSLINC)
    link_directories(SDSLLIB)
    install(DIRECTORY DESTINATION ${OUTPUT})
    set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OUTPUT})
endif()

file(GLOB_RECURSE ProjectFiles *.cc)

set(LIBS roaring sdsl pthread)
set(DEPS CROAR SDSL)
include_directories(${PROJECT_SOURCE_DIR})

add_executable(howdesbt ${ProjectFiles})
if (NOT SDSLLIB OR NOT SDSLINC)
    add_dependencies(howdesbt CROAR SDSL)
else()
    add_dependencies(howdesbt CROAR)
endif()
target_link_libraries(howdesbt ${LIBS})




