#ifndef CALCTHREADFUNCTION_H
#define CALCTHREADFUNCTION_H

#include "common_struct.h"

/**************************
���ƣ�calcThreadFunction.h
�����������̵߳ĺ���ͷ�ļ�
***************************/


/**************************
���ƣ�struct CalcInfo
������������̴߳��ݵĲ����ṹ��
***************************/
struct CalcInfo
{
	int socket;              //�ͻ����׽���
	Element *points;	     //ģ�Ͷ�����Ϣ
	Triangle *triangles;     //ģ��������Ԫ��Ϣ
  Axis_slx *recvPoints;    //��ģ�ͽ��յ���Ϣ
	int points_length;       //���㳤��
	int triangles_length;    //��Ԫ����
	ConfigStruct config;     //�������ò���
};

/**************************
���ƣ�void *recvThreadFunction(void *argv);
�����������̵߳ĺ���
��������Ҫ���̴߳��ݵĲ���
����ֵ����
***************************/
void* calcThreadFunction(void *argv);

#endif
