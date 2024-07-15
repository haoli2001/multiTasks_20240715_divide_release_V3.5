#include "GPUWatchThreadFun.h"

//使用NVIDIA官方提供的nvml库获取GPU状态信息
#include <nvml.h>
#include "common_struct.h"
#include <memory.h>
#include <malloc.h>
#include <string.h>
#include <unistd.h>
#include "socketFunctions.h"

#ifdef __WINDOWS_
#include <Windows.h>
#endif
#ifdef _WIN32
#include <Windows.h>
#endif



void *GPUWatchThreadFun(void *argv)
{
	int socketClient = *(int*)argv;
	nvmlInit();
	unsigned int deviceCount;

	/*@return
		*-\ref NVML_SUCCESS                 if \a deviceCount has been set
		*         -\ref NVML_ERROR_UNINITIALIZED     if the library has not been successfully initialized
		*         -\ref NVML_ERROR_INVALID_ARGUMENT  if \a deviceCount is NULL
		*         -\ref NVML_ERROR_UNKNOWN           on any unexpected error
	*/
	if (NVML_ERROR_UNINITIALIZED == nvmlDeviceGetCount(&deviceCount))
	{
		printf("the library has not been successfully initialized");
		return NULL;
	}
	if (NVML_ERROR_INVALID_ARGUMENT == nvmlDeviceGetCount(&deviceCount))
	{
		printf("\a deviceCount is NULL");
		return NULL;
	}
	if (NVML_ERROR_UNKNOWN == nvmlDeviceGetCount(&deviceCount))
	{
		printf("on any unexpected error");
		return NULL;
	}
	nvmlDevice_t *device = (nvmlDevice_t*)malloc(sizeof(nvmlDevice_t)*deviceCount);

	nvmlUtilization_t *utilization = (nvmlUtilization_t*)malloc(sizeof(nvmlUtilization_t)*deviceCount);

	while (true)
	{
		GPUWatchStruct *gpuwatchstruct = (GPUWatchStruct*)malloc(sizeof(GPUWatchStruct)*deviceCount);
		for (int i = 0; i < deviceCount; i++)
		{
			nvmlDeviceGetHandleByIndex(i, &device[i]);
			nvmlDeviceGetUtilizationRates(device[i], &utilization[i]);
			gpuwatchstruct[i].device_id = i;
			gpuwatchstruct[i].gpu = utilization[i].gpu;
			gpuwatchstruct[i].memory = utilization[i].memory;
			nvmlDeviceGetTemperatureThreshold(device[i], NVML_TEMPERATURE_THRESHOLD_SHUTDOWN, (unsigned int*)&gpuwatchstruct[i].shutdown_temp);
			nvmlDeviceGetTemperatureThreshold(device[i], NVML_TEMPERATURE_THRESHOLD_SLOWDOWN, (unsigned int*)&gpuwatchstruct[i].slowdown_temp);
			nvmlDeviceGetTemperature(device[i], NVML_TEMPERATURE_GPU, (unsigned int*)&gpuwatchstruct[i].temp);
			nvmlMemory_t memory;
			nvmlDeviceGetMemoryInfo(device[i], &memory);
			gpuwatchstruct[i].total = memory.total;
			gpuwatchstruct[i].used = memory.used;
			gpuwatchstruct[i].free = memory.free;
		}
		Frame frame;
		
		for (int i = 0; i < deviceCount; i++)
		{
			strcpy(frame.command, "GPUWatch");
			frame.length = sizeof(GPUWatchStruct);
			memcpy(frame.data, (char*)&gpuwatchstruct[i], sizeof(GPUWatchStruct));
			//send_frame(socketClient, (char*)&frame, sizeof(Frame));
		}
	
        //每隔一秒获取一次信息
#ifdef linux
		sleep(1);
#endif
#ifdef _UNIX
		sleep(1);
#endif
#ifdef __WINDOWS_
		Sleep(1000);
#endif
#ifdef _WIN32
		Sleep(1000);
#endif
	}
}

