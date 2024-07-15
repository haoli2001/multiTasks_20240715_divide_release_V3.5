#include <stdio.h>
#include <pthread.h>
#include <cuda_runtime.h>
#include "helper_cuda.h"
#include <semaphore.h>
#include <pthread.h>
#ifdef linux
#include <sys/types.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <sys/wait.h>
#include <unistd.h>
#include <netinet/tcp.h>

#endif
#ifdef _UNIX
#include <sys/types.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <sys/wait.h>
#include <unistd.h>
#include <netinet/tcp.h>
#endif
#ifdef __WINDOWS_
#include <winsock2.h>

#endif
#ifdef _WIN32
#include <winsock2.h>

#endif


#pragma comment(lib, "Ws2_32.lib")

#include "socketFunctions.h"
#include "common_struct.h"
#include "recvThreadFunction.h"
#include "GPUWatchThreadFun.h"

#define PORT 3490  //定义默认端口为3490端口


//服务器socket和连接的客户端socket 
int socketSevice;


int main()
{
   
	pthread_t recvThread;       //计算线程；
	pthread_t GPUWatchThread;   //GPU 监听线程
	int socketClient[10];       //客户端socket
	int curIndex = 0;           //当前连接的id
    
    
	//socket的初始化，绑定，监听以及等待连接
	init_socket();
    //创建套接字
	socketSevice = create_socket();
    //将创建的套接字绑定到指定端口
	if (-1 == bind_listen(socketSevice, PORT))
	{ 
		printf("bind&listen wrong!\n");
		return 0;
	}
	
    //循环监听连接信息
    while(true) 
    {
        //监听连接，当建立建立后accept_client 函数返回，否则阻塞
        socketClient[curIndex] = accept_client(socketSevice);
        
        int deviceCount;
        //获取服务器GPU设备数
        cudaGetDeviceCount(&deviceCount); 
        for (int dev = 0; dev < deviceCount; dev++)
        {
            cudaSetDevice(dev);
            cudaDeviceProp deviceProp;
            cudaGetDeviceProperties(&deviceProp, dev);
            Frame frame;
            DeviceInfo gpuInfo;
            strcpy(frame.command, "DeviceInfo");
            gpuInfo.deviceID = dev;
            gpuInfo.deviceCount = deviceCount;
            strcpy(gpuInfo.deviceName, deviceProp.name);
            //获取设备每个多处理器的核心数量
            gpuInfo.coresPreMutiprocess = _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor); 
            //多处理器数量 
            gpuInfo.mutiprocessCount = deviceProp.multiProcessorCount;          
            frame.length = sizeof(DeviceInfo);
            memcpy(frame.data, (char*)&gpuInfo, sizeof(gpuInfo));
            //向客户端发送设备信息
            send_frame(socketClient[curIndex], (char*)&frame, sizeof(Frame));     
        }
        pthread_create(&recvThread, NULL, recvThreadFunction, (void*)&(socketClient[curIndex]));    //创建接收数据的线程
        pthread_create(&GPUWatchThread,NULL , GPUWatchThreadFun, (void*)&(socketClient[curIndex]));  //创建GPU实时监控的线程
        curIndex++;
    }
	
    //计算结束，关闭socket
#ifdef linux
	for(int i = curIndex-1;i>=0;i--)
	{	
		shutdown(socketClient[i], SHUT_RDWR);
		close(socketClient[i]);
	}
	shutdown(socketSevice, SHUT_RDWR);
	close(socketSevice);
#endif
#ifdef _UNIX
	for(int i = curIndex-1;i>=0;i--)
	{	
		shutdown(socketClient[i], SHUT_RDWR);
		close(socketClient[i]);
	}
	shutdown(socketClient, SHUT_RDWR);
	shutdown(socketSevice, SHUT_RDWR);
	close(socketClient);
	close(socketSevice);
#endif
#ifdef __WINDOWS_
	for(int i = curIndex-1;i>=0;i--)
	{	
		close(socketClient[i]);
	}
	close(socketSevice);
#endif
#ifdef _WIN32
	for(int i = curIndex-1;i>=0;i--)
	{	
		close(socketClient[i]);
	}
	closesocket(socketSevice);
#endif
	return 0;
}
