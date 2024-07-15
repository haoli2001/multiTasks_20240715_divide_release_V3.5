#ifndef CALCTHREADFUNCTION_H
#define CALCTHREADFUNCTION_H

#include "common_struct.h"

/**************************
名称：calcThreadFunction.h
描述：计算线程的函数头文件
***************************/


/**************************
名称：struct CalcInfo
描述：向计算线程传递的参数结构体
***************************/
struct CalcInfo
{
	int socket;              //客户端套接字
	Element *points;	     //模型顶点信息
	Triangle *triangles;     //模型三角面元信息
  Axis_slx *recvPoints;    //多模型接收点信息
	int points_length;       //顶点长度
	int triangles_length;    //面元长度
	ConfigStruct config;     //解算配置参数
};

/**************************
名称：void *recvThreadFunction(void *argv);
描述：计算线程的函数
参数：需要向线程传递的参数
返回值：无
***************************/
void* calcThreadFunction(void *argv);

#endif
