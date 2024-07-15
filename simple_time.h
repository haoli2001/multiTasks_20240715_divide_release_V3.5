/********************************************************************
	���:       header_file
	����:		2018/10/26
	�ļ���: 	simple_time.h
	����:	    ����һ�δ�����ʱ������һ��simple_time�࣬Windows��Linuxͨ��
*********************************************************************/

#ifdef linux
#include <sys/time.h>
#endif
#ifdef _UNIX
#include <sys/time.h>
#endif
#ifdef __WINDOWS_
#include <time.h>
#endif
#ifdef _WIN32
#include <time.h>
#include <Windows.h>
#endif
#include <stddef.h>
class simple_time
{
public:
	simple_time();
	~simple_time();

	void Time_Start();
	float Time_End();
private:
#ifdef linux
	struct timeval ta, tb;
	double mseca, msecb;
#endif
#ifdef _WIN32
	LARGE_INTEGER  num;
	long long start, end, freq, time1, time2;
#endif
};