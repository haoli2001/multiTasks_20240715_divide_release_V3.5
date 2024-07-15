#ifndef _TREE2VECTOR_H_
#define _TREE2VECTOR_H_
#include "common_struct.h"
/**************************
名称：tree2vector.h
描述：接收数据线程的函数头文件
***************************/



/**************************
名称：KD_Node_V* tree2vector(KD_Node *nodec, int* length);
描述：将树型结构的KDTree转换为连续存储结构的KDTree
参数：KD_Node *nodec:   树状存储结构
      int* length：     转换后的数组长度
返回值：连续存储的KDTree
***************************/

KD_Node_V* tree2vector(KD_Node *nodec, int* length);

#endif