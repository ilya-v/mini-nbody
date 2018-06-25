#ifndef TIMER_H
#define TIMER_H

#include <stdlib.h>
#include <time.h>

static clock_t t0;

void StartTimer()
{
  t0 = clock();
}


double GetTimer()
{
  clock_t t = clock();
  return (t - t0) * 1000.0 / CLOCKS_PER_SEC;
}

#endif // TIMER_H
