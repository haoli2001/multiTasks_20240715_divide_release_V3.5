#ifndef __kd_struct_H__
#define __kd_struct_H__
#define Points_Num 0x130ec//78078//顶点数量
#define Max_Depth 20//最大树高
#define DATA_Type float//坐标信息的精度
#define DataAbs fabs//根据数据e类型调整求绝对值的函数类型
#define Prim_Min 32//节点内最少的面元数量
#define BigNodeNum 256
#define M 3//矩阵维度
#include <stdio.h>
#include <string.h>

#include "common_struct.h"

struct KD_Infomation
{
	long total_node_num;
	long leaf_node_num;
	long empty_node_num;
	long intermediate_node_num;//中间节点
	int depth_max;
	int prim_max;
	int prim_count;
};
void Preprocessing_Triangles(FILE *fp, Triangle *Triangles, Element *Points, int *Node_Num,int *Triangle_Num);
void Build_BigNode(struct KD_Node * treenode, struct Prim_Box *array, struct Prim_Box *out_array, int arr_length, int *out_arr_length);
void Sort_Box(Element points[], Triangle triangle[], Prim_Box *prim, int Triangle_Num);
void Optimization_Rope(KD_Node *treenode);
void KD_Node_init(struct KD_Node *kd_node, Prim_Box *array,int arr_length);
void Destroy_Tree(KD_Node *root);

#endif
