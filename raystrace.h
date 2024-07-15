#ifndef RAYSTRACE_H
#define RAYSTRACE_H

#include "common_struct.h"

/**************************
���ƣ�raystrace.h
����������׷�ٺ���ͷ�ļ�
***************************/




/**************************
���ƣ�void allraystrace_v2(Direction *d_rays, Square *d_squares, int width, int height, KD_Node_V *d_root, Prim_Box *d_array, Element *d_points, Triangle *d_triangles,
	int* d_DivRayTubeNum, int* DivRayTubeNum, int* d_sum_gmem, int* d_sum_Gmem, int* d_squares_pred, Axis direction, float angle);
����������׷�ٺ���
������Direction *d_rays����������
      Square *d_squares���׾�����
      int width, int height��   �׾�����
      KD_Node_V *d_root��       ���Դ洢KDTree
      Prim_Box *d_array, Element *d_points, Triangle *d_triangles�� ��Χ�����ݣ�ģ�ͽڵ����ݣ�ģ������������
      int* d_DivRayTubeNum, int* DivRayTubeNum�� 
      int* d_sum_gmem, int* d_sum_Gmem ��
      int* d_squares_pred��
      Axis direction, float angle��
      
����ֵ����
***************************/
void allraystrace_v2(Direction *d_rays, Square *d_squares, int width, int height, KD_Node_V *d_root, Prim_Box *d_array, Element *d_points, Triangle *d_triangles,
	int* d_DivRayTubeNum, int* DivRayTubeNum, int* d_sum_gmem, int* d_sum_Gmem, int* d_squares_pred, Axis direction, float angle, float water_line);

void allraystrace_DivRayFirst(Direction *d_rays2, Square *d_squares2, int d_totalDivRayNum, KD_Node_V *d_root, Prim_Box *d_array, Element *d_points, Triangle *d_triangles,
	int* d_DivRayTubeNum, int* DivRayTubeNum, int* d_sum_gmem, int* d_sum_Gmem, int* d_squares_pred, Axis direction, float angle, float water_line);

void allraystrace_DivRaySecond(Direction *d_rays2, Square *d_squares2, int d_totalDivRayNum, KD_Node_V *d_root, Prim_Box *d_array, Element *d_points, Triangle *d_triangles,
	int* d_DivRayTubeNum, int* DivRayTubeNum, int* d_sum_gmem, int* d_sum_Gmem, int* d_squares_pred, Axis direction, float angle, float water_line);

void allraystrace_DivRayThird(Direction *d_rays2, Square *d_squares2, int d_totalDivRayNum, KD_Node_V *d_root, Prim_Box *d_array,
	Element *d_points, Triangle *d_triangles, Axis direction, float water_line);

    
#endif