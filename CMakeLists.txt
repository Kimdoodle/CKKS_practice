# CMakeList.txt: SEAL_VS에 대한 CMake 프로젝트, 여기에 소스를 포함하고
# 프로젝트 특정 논리를 정의합니다.
#

# 이 프로젝트의 실행 파일에 소스를 추가합니다.
add_executable (SEAL_VS "SEAL_VS.cpp" "SEAL_VS.h" "compare.cpp" "compare.h" "arithmetic_seal.cpp" "arithmetic_seal.h" "polymath.cpp" "polymath.h" "print.cpp" "print.h" "function_seal.cpp" "function_seal.h")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET SEAL_VS PROPERTY CXX_STANDARD 20)
endif()

# TODO: 필요한 경우 테스트를 추가하고 대상을 설치합니다.

set(SEAL_DIR "C:/Program Files (x86)/SEAL")

# SEAL 라이브러리 찾기
find_package(SEAL REQUIRED PATHS "${SEAL_DIR}" NO_DEFAULT_PATH)

# SEAL 라이브러리 링크
target_link_libraries(SEAL_VS PRIVATE SEAL::seal)

# Debug 모드에서 _ITERATOR_DEBUG_LEVEL=2 설정
target_compile_definitions(SEAL_VS PRIVATE "_ITERATOR_DEBUG_LEVEL=2")

# Debug 모드에서 Multi-threaded Debug (/MDd) 강제 설정
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(SEAL_VS PRIVATE /MDd)
endif()