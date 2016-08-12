CMAKE_MINIMUM_REQUIRED(VERSION 3.6)

PROJECT("2-3D phase congruency")

SET(INCLUDE_OUTPUT_PATH
    ${CMAKE_BINARY_DIR}/include
    CACHE PATH "Single output directory for all includes.")

SET(EXECUTABLE_OUTPUT_PATH
    ${CMAKE_BINARY_DIR}/bin
    CACHE PATH "Single output directory for building all executables.")

SET(LIBRARY_OUTPUT_PATH
    ${CMAKE_BINARY_DIR}/lib
    CACHE PATH "Single output directory for building all libraries.")

INCLUDE_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party
    ${CMAKE_CURRENT_SOURCE_DIR}/utils
    ${CMAKE_CURRENT_SOURCE_DIR}/filter_bank
    ${CMAKE_CURRENT_SOURCE_DIR}/phase_congruency)

ADD_SUBDIRECTORY(utils)
ADD_SUBDIRECTORY(filter_bank)
ADD_SUBDIRECTORY(phase_congruency)
