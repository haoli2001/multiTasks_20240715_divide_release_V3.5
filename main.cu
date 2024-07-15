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

#define PORT 3490  //����Ĭ�϶˿�Ϊ3490�˿�


//������socket�����ӵĿͻ���socket 
int socketSevice;


int main()
{
   
	pthread_t recvThread;       //�����̣߳�
	pthread_t GPUWatchThread;   //GPU �����߳�
	int socketClient[10];       //�ͻ���socket
	int curIndex = 0;           //��ǰ���ӵ�id
    
    
	//socket�ĳ�ʼ�����󶨣������Լ��ȴ�����
	init_socket();
    //�����׽���
	socketSevice = create_socket();
    //���������׽��ְ󶨵�ָ���˿�
	if (-1 == bind_listen(socketSevice, PORT))
	{ 
		printf("bind&listen wrong!\n");
		return 0;
	}
	
    //ѭ������������Ϣ
    while(true) 
    {
        //�������ӣ�������������accept_client �������أ���������
        socketClient[curIndex] = accept_client(socketSevice);
        
        int deviceCount;
        //��ȡ������GPU�豸��
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
            //��ȡ�豸ÿ���ദ�����ĺ�������
            gpuInfo.coresPreMutiprocess = _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor); 
            //�ദ�������� 
            gpuInfo.mutiprocessCount = deviceProp.multiProcessorCount;          
            frame.length = sizeof(DeviceInfo);
            memcpy(frame.data, (char*)&gpuInfo, sizeof(gpuInfo));
            //��ͻ��˷����豸��Ϣ
            send_frame(socketClient[curIndex], (char*)&frame, sizeof(Frame));     
        }
        pthread_create(&recvThread, NULL, recvThreadFunction, (void*)&(socketClient[curIndex]));    //�����������ݵ��߳�
        pthread_create(&GPUWatchThread,NULL , GPUWatchThreadFun, (void*)&(socketClient[curIndex]));  //����GPUʵʱ��ص��߳�
        curIndex++;
    }
	
    //����������ر�socket
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
