#include <time.h>

double get_timediff_mus(struct timespec start_time, struct timespec stop_time) {
  return difftime(stop_time.tv_sec, start_time.tv_sec) * 1000000 +
         (stop_time.tv_nsec - start_time.tv_nsec) / 1000.0;
}
