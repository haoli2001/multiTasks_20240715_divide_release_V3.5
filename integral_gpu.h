#ifndef __INTEGRAL_GPU_H__
#define __INTEGRAL_GPU_H__
#include "common_struct.h"

typedef struct vector
{
	float x;
	float y;
	float z;
}Vector;//12Bytes

typedef struct fushu
{
	float re;
	float im;
}comp;

typedef struct ouput
{
	float re;
	float im;
	int triangle_index;
}ReimOutput;

//typedef struct single_ray
//{
//	//float p[4];					//交点
//	float distance;				//交点与发射点距离
//}SingleRayInfo;

typedef struct point_coordinate
{
	float p[3];
}PointCoor3;//12Bytes

typedef struct point_coordinates
{
	float p[4];
}PointCoor;//16Bytes

typedef struct matrix_struct
{
	float p[12];
}MatStruct;//48Bytes

//单射线保存的信息
typedef struct beam_ray
{
	PointCoor launch_point;
	PointCoor ray_index[4];		//射线对应的4个顶点
	PointCoor point_2D[5];
	float p_cent_distance;
	PointCoor3 parameter;
	float gorden;
	float Z0;
 
  float gorden_re;
  float gorden_im;
	
}RayBeamInfo;//176Bytes

 
 comp sound_field_integral_gpu(Direction *d_rays, Square *d_squares, float wavelength, int raysBeamNum, RayBeamInfo* d_effrays, Vector* d_center, Vector* d_axis, 
	MatStruct* d_transMat, ReimOutput* d_reim, float* d_sum_re, float* d_sum_im, float fi, float si,  ConfigStruct* config);//20210308 wangying 20210831姬梓遇
 

comp sound_field_integral_gpu_DivRay(Direction *d_rays, Square *d_squares, float wavelength, int raysBeamNum,
	RayBeamInfo* d_effrays, Vector* d_center, Vector* d_axis, MatStruct* d_transMat, ReimOutput* d_reim, float* d_sum_re, float* d_sum_im, float fi, float si,  ConfigStruct* config);

comp sound_field_integral_gpu_DivRay3(Direction *d_rays, Square *d_squares, float wavelength, int raysBeamNum,
	RayBeamInfo* d_effrays, Vector* d_center, Vector* d_axis, MatStruct* d_transMat, ReimOutput* d_reim, float* d_sum_re, float* d_sum_im, float fi, float si,  ConfigStruct* config);

float TS_compute(float far_dis,comp sum, float wavelength);

void DivRayTube(Direction *d_rays, Square *d_squares, Direction *d_rays2, Square *d_squares2, int* d_DivRayTubeNum, int DivRayTubeNum, int* d_sum_Gmem,
	int* d_squares_pred, int RayTubeNum, float lmd, float fi, Axis direction);

#endif
