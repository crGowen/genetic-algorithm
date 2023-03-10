cmake_minimum_required(VERSION 3.16)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED True)

set(OUTDIR "${CMAKE_SOURCE_DIR}/output")

# Print building info
message("\nConfiguration:")
message("\tOS: ${CMAKE_HOST_SYSTEM_NAME}")
message("\tPlatform: ${CMAKE_GENERATOR_PLATFORM}")
message("\tDirectory: ${CMAKE_SOURCE_DIR}")
message("\tProcessor: ${CMAKE_HOST_SYSTEM_PROCESSOR}")

# Building genalg library
project(genalg)
add_library(genalg SHARED ./genalg/main.c)
target_compile_options(genalg PRIVATE -Wall -O2)

set(GENALG_LIB_INCLUDE "${CMAKE_SOURCE_DIR}/genalg/include")
target_include_directories(genalg PRIVATE "${GENALG_LIB_INCLUDE}")
target_compile_definitions(genalg PRIVATE GENALG_EXPORT)

add_custom_command(
    TARGET genalg POST_BUILD
    COMMAND "${CMAKE_COMMAND}" -E copy
        "$<TARGET_LINKER_FILE:genalg>"
        "${OUTDIR}/lib/$<TARGET_LINKER_FILE_NAME:genalg>"
    COMMAND "${CMAKE_COMMAND}" -E copy_directory
        "${GENALG_LIB_INCLUDE}"
        "${OUTDIR}/include"
)
if (${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Windows")
    # .dll copy for Windows
    add_custom_command(
        TARGET genalg POST_BUILD
        COMMAND "${CMAKE_COMMAND}" -E copy
            "$<TARGET_FILE:genalg>"
            "${OUTDIR}/$<TARGET_FILE_NAME:genalg>"
    )
endif()

# Building usage example
project(usage_example)
add_executable(usage_example ./usage_example/main.c)
target_compile_options(usage_example PRIVATE -Wall -O2)

add_dependencies(usage_example genalg)
target_include_directories(usage_example PRIVATE "${OUTDIR}/include")
target_link_libraries(usage_example PRIVATE "${OUTDIR}/lib/$<TARGET_LINKER_FILE_NAME:genalg>")
add_custom_command(
    TARGET usage_example POST_BUILD
    COMMAND "${CMAKE_COMMAND}" -E copy
        "$<TARGET_FILE:usage_example>"
        "${OUTDIR}/$<TARGET_FILE_NAME:usage_example>"
)

# Building usage example
project(test)
add_executable(test ./test/main.c)
target_compile_options(test PRIVATE -Wall -O2)

add_dependencies(test genalg)
target_include_directories(test PRIVATE "${OUTDIR}/include")
target_link_libraries(test PRIVATE "${OUTDIR}/lib/$<TARGET_LINKER_FILE_NAME:genalg>")
add_custom_command(
    TARGET test POST_BUILD
    COMMAND "${CMAKE_COMMAND}" -E copy
        "$<TARGET_FILE:test>"
        "${OUTDIR}/$<TARGET_FILE_NAME:test>"
)