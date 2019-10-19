#include "Timer.h"

#include <windows.h>

Timer::Timer()
{
  LARGE_INTEGER freq;
  QueryPerformanceFrequency(&freq);
  m_frequency = (double)freq.QuadPart;
}

double Timer::getTimeStamp()
{
  LARGE_INTEGER query_ticks;
  QueryPerformanceCounter(&query_ticks);
  return query_ticks.QuadPart/m_frequency;
}
