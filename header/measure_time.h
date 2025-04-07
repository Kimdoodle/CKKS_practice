#ifndef MEASURE_TIME_H
#define MEASURE_TIME_H

#include "SEAL_VS.h"

time_point<high_resolution_clock> cur_time();
void calculate_time(time_point<high_resolution_clock> time1, time_point<high_resolution_clock> time2);

#endif //MEASURE_TIME_H