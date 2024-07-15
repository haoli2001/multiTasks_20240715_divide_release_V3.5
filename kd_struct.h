#ifndef __kd_struct_H__
#define __kd_struct_H__
#define Points_Num 0x130ec//78078//��������
#define Max_Depth 20//�������
#define DATA_Type float//������Ϣ�ľ���
#define DataAbs fabs//��������e���͵��������ֵ�ĺ�������
#define Prim_Min 32//�ڵ������ٵ���Ԫ����
#define BigNodeNum 256
#define M 3//����ά��
#include <stdio.h>
#include <string.h>

#include "common_struct.h"

struct KD_Infomation
{
	long total_node_num;
	long leaf_node_num;
	long empty_node_num;
	long intermediate_node_num;//�м�ڵ�
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
