// SEAL_VS.h: 표준 시스템 포함 파일
// 또는 프로젝트 특정 포함 파일이 들어 있는 포함 파일입니다.

#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <random>
#include <memory>
#include <stdexcept>
#include <iomanip>
#include <chrono>
#include "seal/seal.h"

using namespace std;
using namespace seal;
using namespace std::chrono;

// TODO: 여기서 프로그램에 필요한 추가 헤더를 참조합니다.
#include "CKKS_class.h"
#include "polymath.h"
#include "function_plain.h"
#include "function_seal.h"
#include "print.h"
#include "measure_time.h"

