cmake_minimum_required (VERSION 3.8)

if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

set(PROJECT_NAME "SEAL_SGN")
    project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB SOURCES "*.cpp")
add_executable(${PROJECT_NAME} ${SOURCES})

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET SEAL_SGN PROPERTY CXX_STANDARD 20)
endif()


if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(SEAL_LIB_DIR "${CMAKE_SOURCE_DIR}/SEAL/build/lib/Debug")
else()
    set(SEAL_LIB_DIR "${CMAKE_SOURCE_DIR}/SEAL/build/lib/Release")
endif()

set(SEAL_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/SEAL/build/native/src")

# include & lib path 
target_include_directories(${PROJECT_NAME} PRIVATE
    "${CMAKE_SOURCE_DIR}/SEAL/native/src"
)
target_include_directories(${PROJECT_NAME} PRIVATE
    "${CMAKE_SOURCE_DIR}/SEAL/build/native/src"
)
target_include_directories(${PROJECT_NAME} PRIVATE
    "${CMAKE_SOURCE_DIR}/SEAL/build/thirdparty//msgsl-src/include"
)

target_link_directories(${PROJECT_NAME} PRIVATE ${SEAL_LIB_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/header)
target_link_libraries(${PROJECT_NAME} PRIVATE seal-4.1)