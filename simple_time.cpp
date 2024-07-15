#include "simple_time.h"


simple_time::simple_time()
{

}

simple_time::~simple_time()
{

}

void simple_time::Time_Start()
{
#ifdef linux
	gettimeofday(&ta, NULL);
	mseca = ta.tv_sec*1000.0 + ta.tv_usec / 1000.0;
#endif
#ifdef _WIN32
	QueryPerformanceFrequency(&num);
	freq = num.QuadPart;
	QueryPerformanceCounter(&num);
	start = num.QuadPart;
#endif
}

float simple_time::Time_End()
{
#ifdef linux
	gettimeofday(&tb, NULL);
	msecb = tb.tv_sec*1000.0 + tb.tv_usec / 1000.0;
	msecb -= mseca;
	return msecb;
#endif
#ifdef _WIN32
	QueryPerformanceCounter(&num);
	end = num.QuadPart;
	time1 = (end - start) * 1000.0 / freq;
	return (end - start) * 1000.0 / freq;
#endif
}