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

if (NOT SDSLLIB OR NOT SDSLINC OR NOT KMTRICKS_INC OR NOT KMTRICKS_LIB)
    message("-DSDSLLIB=path -DSDSLINC=path -DKMTRICKS_INC=path -DKMTRICKS_LIB=path are required")
else()

    if(NOT OUTPUT)
        set(OUTPUT ${PROJECT_SOURCE_DIR}/bin)
    endif()

    include_directories(${SDSLINC} ${KMTRICKS_INC})
    link_directories(${SDSLLIB} ${KMTRICKS_LIB})

    install(DIRECTORY DESTINATION ${OUTPUT})
    set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OUTPUT})

    file(GLOB_RECURSE ProjectFiles *.cc)

    set(LIBS kmtricks roaring sdsl pthread z)
    set(DEPS CROAR SDSL)
    include_directories(${PROJECT_SOURCE_DIR})

    set(CMAKE_INSTALL_RPATH ${PROJECT_BINARY_DIR}/CROAR/src/CROAR-build/)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
    add_executable(howdesbt ${ProjectFiles})
    add_dependencies(howdesbt CROAR)
    target_link_libraries(howdesbt ${LIBS})

endif()




