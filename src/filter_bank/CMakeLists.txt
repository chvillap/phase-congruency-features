CMAKE_MINIMUM_REQUIRED(VERSION 3.6)

PROJECT("filter_bank")

# -------------------------------------------------------------------------------------------------

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

# -------------------------------------------------------------------------------------------------

SET(LIBRARY_NAME
    "bip_filter_bank"
)
SET(LIBRARY_SOURCE_FILES
    "log_gabor_filter_bank.cpp"
)
ADD_LIBRARY(
    ${LIBRARY_NAME}
    ${LIBRARY_SOURCE_FILES}
)
TARGET_LINK_LIBRARIES(
    ${LIBRARY_NAME}
    ${ITK_LIBRARIES}
    bip_utils
)
TARGET_INCLUDE_DIRECTORIES(
    ${LIBRARY_NAME}
    PUBLIC
    ${CMAKE_SOURCE_DIR}/${PROJECT_NAME}
)
SET_PROPERTY(
    TARGET ${LIBRARY_NAME}
    PROPERTY CXX_STANDARD 11)
SET_PROPERTY(
    TARGET ${LIBRARY_NAME}
    PROPERTY CXX_STANDARD_REQUIRED ON)

# -------------------------------------------------------------------------------------------------

SET(PROGRAM_NAME
    "log_gabor_filter_bank"
)
SET(PROGRAM_SOURCE_FILES
    "main.cpp"
)
ADD_EXECUTABLE(
    ${PROGRAM_NAME}
    ${PROGRAM_SOURCE_FILES}
)
TARGET_LINK_LIBRARIES(
    ${PROGRAM_NAME}
    ${ITK_LIBRARIES}
    bip_utils
    ${LIBRARY_NAME}
)
SET_PROPERTY(
    TARGET ${PROGRAM_NAME}
    PROPERTY CXX_STANDARD 11)
SET_PROPERTY(
    TARGET ${PROGRAM_NAME}
    PROPERTY CXX_STANDARD_REQUIRED ON)
