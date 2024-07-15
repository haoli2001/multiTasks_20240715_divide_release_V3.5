#ifndef _VIRTUALFACE_CUDA_H
#define _VIRTUALFACE_CUDA_H

#include "common_struct.h"

struct Radius						//������������
{
	float Xr;
	float Yst;
	float Zfi;
};

struct DynamicPlane
{
	float st_min;
	float fi_max;
	float ratio;
	int width;
	int height;
	float GPUnum;
	bool flag;
};

struct BinaryTimeTree
{
	bool flag;
	int GPUnum;
	float runtime;
	float ratio;
	float AveTimeLeft;
	float AveTimeRight;
	BinaryTimeTree *leftchild;
	BinaryTimeTree *rightchild;
};

Axis dSphericaltoRectangular(Radius a);
void getWidthHeight(float far_dis,Element *points, int num, int *width_calc, int *height_calc,float st, float fi, float lmd_calc, int *dheight_calc, float *de_st_min, float *de_fi_max, int device_num);
void create_virtualface_gpu(Direction *rays2, Square *squares2, int dwidth, int dheight, float lmd, Radius direction_radius,
	float e_st_min, float e_fi_max, int device_id);
void free_virtualface(Direction *rays, Square *squares);
void getWidthHeight_Sun(float far_dis,Box b, int *width, int *height, int st, float fi, float lmd, int *dheight, float *de_st_min, float *de_fi_max, int device_num);
void copy2host(Direction *d_rays, Square *d_squares, Direction **rays, Square **squares, int width, int height);
//DynamicPlane* ConstructVirtualFace(DynamicPlane* dData, BinaryTimeTree** preangletime, int DeviceCount, float e_st_min, float e_fi_max, int width, int height, float lmd);
BinaryTimeTree** ConstructTimeTree(DynamicPlane* tmp, float *runtime, int DeviceCount);
void ConstructVirtualFace(DynamicPlane* array, DynamicPlane* dData, BinaryTimeTree** pre_angle_time, int DeviceCount, float e_st_min, float e_fi_max, int width, int height, float lmd);
#endif
