#include "recvThreadFunction.h"

#include <stdio.h>
#include <semaphore.h>
#include <pthread.h>
#include<signal.h>
#include <unistd.h>
#include <string.h>
#include <malloc.h>
#include "calcThreadFunction.h"
#include "common_struct.h"
#include "socketFunctions.h"
extern pthread_mutex_t socket_mutex;//lihao 20240711传输锁
void *recvThreadFunction(void *argv)
{
	CalcInfo calcInfo;                //解算参数配置
	int socketClient = *(int*)argv;   //接收套接字
    
	calcInfo.socket = socketClient; 

	pthread_t calcThread;             //计算线程    
	memset(&calcThread, 0, sizeof(pthread_t));
    
    //循环接收数据
	while (true)
	{
		Frame frame;
        //读取一个数据帧
		recv_data(socketClient, (char*)&frame, sizeof(Frame));

		if (!strcmp(frame.command, "Triangles"))
		{

			//接收triangles数据
			if(calcInfo.triangles != NULL)
			{
				free(calcInfo.triangles);
			}
			printf("triangles data recv start!\n");
			calcInfo.triangles_length = frame.length / sizeof(Triangle);
			calcInfo.triangles = (Triangle *)malloc(frame.length);
			recv_data(socketClient, (char*)calcInfo.triangles, frame.length);
			printf("triangles data recv over!\n");
		}
		if (!strcmp(frame.command, "Elements"))
		{
			//接收Element数据
			if(calcInfo.points != NULL)
			{
				free(calcInfo.points);
			}
			printf("elements data recv start!\n");
			calcInfo.points_length = frame.length / sizeof(Element);
			calcInfo.points = (Element *)malloc(frame.length);
			recv_data(socketClient, (char*)calcInfo.points, frame.length);
			printf("elements data recv over!\n");

			//数据接收完毕后，发送等待配置帧
			Frame frame;
			strcpy(frame.command, "WaitForConfig");
			frame.length = 0;
#ifdef linux
			pthread_mutex_lock(&socket_mutex);//lihao 20240711 
#endif
			send_frame(socketClient, (char*)&frame, sizeof(Frame));
#ifdef linux
			pthread_mutex_unlock(&socket_mutex);//lihao 20240711  发送数据时一直占有锁
#endif
		}
		if (!strcmp(frame.command, "RecvPoints"))
		{

			//接收RecvPoints数据
			if(calcInfo.recvPoints != NULL)
			{
				free(calcInfo.recvPoints);
			}
			printf("recvPoints data recv start!\n");
			calcInfo.recvPoints = (Axis_slx *)malloc(frame.length);
			recv_data(socketClient, (char*)calcInfo.recvPoints, frame.length);
			printf("recvPoints data recv over!\n");
		}
		if (!strcmp(frame.command, "Configuration"))
		{
			//配置并开始
			memcpy(&calcInfo.config, frame.data, sizeof(ConfigStruct));

			//开始执行计算线程，当计算线程正在执行时，则先关闭线程后再重新执行
			if (calcThread!=0 && pthread_kill(calcThread, 0) == 0)
			{
				pthread_cancel(calcThread);
				pthread_join(calcThread,NULL);
				printf("restart");
			}
			pthread_create(&calcThread, NULL, calcThreadFunction, (void*)&calcInfo);
		}
		if (!strcmp(frame.command, "Stop"))
		{
            //停止命令，退出线程
			if (calcThread!=0 && pthread_kill(calcThread, 0) == 0)
			{
				pthread_cancel(calcThread);
				pthread_join(calcThread,NULL);
			}
			//2020.3.24 jzy
			memset(&calcThread, 0, sizeof(pthread_t));
		}
		if (!strcmp(frame.command, "Exit"))
		{
            //断开连接，退出接收循环
			if (calcThread!=0 && pthread_kill(calcThread, 0) == 0)
			{
				pthread_cancel(calcThread);
				pthread_join(calcThread,NULL);
			}
			break;
		}
	}
	return NULL;
}
