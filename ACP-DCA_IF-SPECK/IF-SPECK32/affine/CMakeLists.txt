CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(ACPDCA)

include_directories(include)

aux_source_directory(src/m4ri ACPDCA_SRC)

add_library(ACPDCA_LIB ${ACPDCA_SRC})

add_executable(ACPDCA test/white_box_arx.c)

target_link_libraries(ACPDCA ACPDCA_LIB)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -O2")