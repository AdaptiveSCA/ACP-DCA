CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(ACPDCA)

include_directories(include)

aux_source_directory(src/WBMatrix ACPDCA_SRC)
aux_source_directory(src ACPDCA_SRC)

add_library(ACPDCA_LIB ${ACPDCA_SRC})

add_executable(ACPDCA test/attack.c)

target_link_libraries(ACPDCA ACPDCA_LIB m)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -O2")