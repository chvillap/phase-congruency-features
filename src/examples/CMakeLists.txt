CMAKE_MINIMUM_REQUIRED(VERSION 3.6)

PROJECT("examples")

# -------------------------------------------------------------------------------------------------

SET(PROGRAM_NAME
    "example1"
)
SET(PROGRAM_SOURCE_FILES
    "example1.cpp"
)
ADD_EXECUTABLE(
    ${PROGRAM_NAME}
    ${PROGRAM_SOURCE_FILES}
)
TARGET_LINK_LIBRARIES(
    ${PROGRAM_NAME}
    bip_utils
    bip_filter_bank
)
SET_PROPERTY(
    TARGET ${PROGRAM_NAME}
    PROPERTY CXX_STANDARD 14)
SET_PROPERTY(
    TARGET ${PROGRAM_NAME}
    PROPERTY CXX_STANDARD_REQUIRED ON)
