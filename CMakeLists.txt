cmake_minimum_required(VERSION 3.7)

project(Laplace)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Werror -fdump-rtl-expand -lm")
SET( CMAKE_EXPORT_COMPILE_COMMANDS ON )

set(SRCS
    main.c
    functions.c
    )

set(HEADERS
    defs.h
    )

IF( EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json" )
    EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
        ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
        )
ENDIF()

add_executable(main ${SRCS} ${HEADERS})
TARGET_LINK_LIBRARIES(main m)
