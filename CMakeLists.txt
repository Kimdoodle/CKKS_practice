# CMakeList.txt: SEAL_SGN_TEMP에 대한 CMake 프로젝트, 여기에 소스를 포함하고
# 프로젝트 특정 논리를 정의합니다.
#
cmake_minimum_required (VERSION 3.8)

# 지원되는 경우 MSVC 컴파일러에 대해 핫 다시 로드 사용하도록 설정합니다.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

set(PROJECT_NAME "SEAL_SGN_TEMP")
    project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 이 프로젝트의 실행 파일에 소스를 추가합니다.
file(GLOB SOURCES "*.cpp")
add_executable(${PROJECT_NAME} ${SOURCES})

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET SEAL_SGN_TEMP PROPERTY CXX_STANDARD 20)
endif()

# TODO: 필요한 경우 테스트를 추가하고 대상을 설치합니다.
# === SEAL 라이브러리 설정 ===
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(SEAL_LIB_DIR "${CMAKE_SOURCE_DIR}/SEAL/build/lib/Debug")
else()
    set(SEAL_LIB_DIR "${CMAKE_SOURCE_DIR}/SEAL/build/lib/Release")
endif()

set(SEAL_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/SEAL/build/native/src")

# include & lib path 추가
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


# link CPET_LIB library.
target_link_libraries(${PROJECT_NAME} PRIVATE seal-4.1)