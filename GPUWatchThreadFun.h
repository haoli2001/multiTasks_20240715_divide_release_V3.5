#ifndef GPUWATCHTHREADFUN_H
#define GPUWATCHTHREADFUN_H


/**************************
名称：GPUWatchThreadFun.h
描述：GPU监控线程的函数头文件，使用NVIDIA官方提供的nvml库获取GPU状态信息
***************************/



/**************************
名称：void *GPUWatchThreadFun(void *argv);
描述：GPU监控线程的函数
参数：监控线程传递的参数
返回值：无
***************************/
void *GPUWatchThreadFun(void *argv);
#endif
