#include "header/SEAL_MM.h"

//return current time.
time_point<high_resolution_clock> cur_time()
{
	return high_resolution_clock::now();
}

//calculate difference between two times.
void calculate_time(time_point<high_resolution_clock> time1, time_point<high_resolution_clock> time2)
{
    auto duration_ms = duration_cast<milliseconds>(time2 - time1);

    if (duration_ms.count() >= 1000) { 
        auto duration_s = duration_cast<duration<double>>(time2 - time1);
        cout << "Execution time: " << duration_s.count() << " s" << endl;
    }
    else {
        cout << "Execution time: " << duration_ms.count() << " ms" << endl;
    }
}




