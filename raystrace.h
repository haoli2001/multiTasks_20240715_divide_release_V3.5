#ifndef RAYSTRACE_H
#define RAYSTRACE_H

#include "common_struct.h"

/**************************
名称：raystrace.h
描述：射线追踪函数头文件
***************************/




/**************************
名称：void allraystrace_v2(Direction *d_rays, Square *d_squares, int width, int height, KD_Node_V *d_root, Prim_Box *d_array, Element *d_points, Triangle *d_triangles,
	int* d_DivRayTubeNum, int* DivRayTubeNum, int* d_sum_gmem, int* d_sum_Gmem, int* d_squares_pred, Axis direction, float angle);
描述：射线追踪函数
参数：Direction *d_rays：射线数据
      Square *d_squares：孔径数据
      int width, int height：   孔径面宽高
      KD_Node_V *d_root：       线性存储KDTree
      Prim_Box *d_array, Element *d_points, Triangle *d_triangles： 包围盒数据，模型节点数据，模型三角形数据
      int* d_DivRayTubeNum, int* DivRayTubeNum： 
      int* d_sum_gmem, int* d_sum_Gmem ：
      int* d_squares_pred：
      Axis direction, float angle：
      
返回值：无
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