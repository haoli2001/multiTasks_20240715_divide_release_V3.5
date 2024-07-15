#include "virtualface_gpu.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "handlerror.h"
#include "device_functions.h"

#define TILE_DIM 16

Radius AxistoRadius_gpu(Axis a, Radius b)		//������ֱ������ϵ���������ת��
{
	Radius c;
	double st, fi;					//st=��,fi=��
	double pi = acos(-1.0);
	st = b.Yst * pi / 180;
	fi = b.Zfi * pi / 180;
	c.Xr = sin(st)*cos(fi)*a.x + sin(st)*sin(fi)*a.y + cos(st)*a.z;
	c.Yst = cos(st)*cos(fi)*a.x + cos(st)*sin(fi)*a.y - sin(st)*a.z;
	c.Zfi = -sin(fi)*a.x + cos(fi)*a.y;
	return c;
}

Axis dSphericaltoRectangular(Radius a)				//���꣬�����굽ֱ������ϵ��ת��
{
	double pi = acos(-1.0);
	Axis b;
	double st, fi;					//st=��,fi=��
	st = a.Yst * pi / 180;
	fi = a.Zfi * pi / 180;

	b.x = a.Xr*sin(st)*cos(fi);
	b.y = a.Xr*sin(st)*sin(fi);
	b.z = a.Xr*cos(st);

	return b;
}

__device__ Axis dSphericaltoRectangular_gpu(Radius a)				//���꣬�����굽ֱ������ϵ��ת��
{
	double pi = acos(-1.0);
	Axis b;
	double st, fi;					//st=��,fi=��
	st = a.Yst * pi / 180;
	fi = a.Zfi * pi / 180;

	b.x = a.Xr*sin(st)*cos(fi);
	b.y = a.Xr*sin(st)*sin(fi);
	b.z = a.Xr*cos(st);

	return b;
}

__global__ void dAxistoRadius_gpu(Radius* d_spherical, Element* d_point, Radius b, int num)		//������ֱ������ϵ���������ת��
{
	double pi = acos(-1.0);
	double st, fi;					//st=��,fi=��
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	int tid = threadIdx.x;
	if (idx < num)
	{
		st = b.Yst * pi / 180;
		fi = b.Zfi * pi / 180;

		d_spherical[idx].Xr = sin(st)*cos(fi)*d_point[idx].point[0] + sin(st)*sin(fi)*d_point[idx].point[1] + cos(st)*d_point[idx].point[2];
		d_spherical[idx].Yst = cos(st)*cos(fi)*d_point[idx].point[0] + cos(st)*sin(fi)*d_point[idx].point[1] - sin(st)*d_point[idx].point[2];
		d_spherical[idx].Zfi = -sin(fi)*d_point[idx].point[0] + cos(fi)*d_point[idx].point[1];
	}

}


__device__ Axis dRadiustoAxis_gpu(Radius a, Radius b)		//�����������굽ֱ������ϵ��ת��
{
	double pi = acos(-1.0);
	Axis c;
	double st, fi;					//st=��,fi=��
	st = b.Yst * pi / 180;
	fi = b.Zfi * pi / 180;
	c.x = sin(st)*cos(fi)*a.Xr + cos(st)*cos(fi)*a.Yst - sin(fi)*a.Zfi;
	c.y = sin(st)*sin(fi)*a.Xr + cos(st)*sin(fi)*a.Yst + cos(fi)*a.Zfi;
	c.z = cos(st)*a.Xr - sin(st)*a.Yst;
	return c;
}

__device__ Axis_slx dRadiustoAxis_gpu_slx(Radius a, Radius b)		//�����������굽ֱ������ϵ��ת��
{
	double pi = acos(-1.0);
	Axis_slx c;
	double st, fi;					//st=��,fi=��
	st = b.Yst * pi / 180;
	fi = b.Zfi * pi / 180;
	c.p[0] = sin(st)*cos(fi)*a.Xr + cos(st)*cos(fi)*a.Yst - sin(fi)*a.Zfi;
	c.p[1] = sin(st)*sin(fi)*a.Xr + cos(st)*sin(fi)*a.Yst + cos(fi)*a.Zfi;
	c.p[2] = cos(st)*a.Xr - sin(st)*a.Yst;
	return c;
}

__global__ void dcreate_virtualface_gpu(Radius direction_radius, Direction *rays, Square *squares, int height, int width, float e_st_min, float e_fi_max, float lmd)
{
	Radius p_cent, p_corner;

	//direction_p = dSphericaltoRectangular_gpu(direction);//等相位面的法向量
	p_cent.Xr = p_corner.Xr = direction_radius.Xr;

	int idx = blockIdx.x * TILE_DIM + threadIdx.x;
	int idy = blockIdx.y * TILE_DIM + threadIdx.y;
	int id = idy * width + idx;
	//每根声线管束的五个射线编号
	if (idx < width && idy < height)
	{
		squares[id].CornerRayIndex.x = idy * width + idy + idx + 1;
		squares[id].CornerRayIndex.y = idy * width + idy + idx;
		squares[id].CornerRayIndex.z = (idy + 1) * width + idy + idx + 1;
		squares[id].CornerRayIndex.w = (idy + 1) * width + idy + idx + 2;
		squares[id].CenterRayIndex = (width + 1) * (height + 1) + idy * width + idx;
		squares[id].right = false;
	}

	if (idx < width && idy < height)//声线管束的中心射线
	{
		p_cent.Yst = e_st_min + (idy + 0.5)*lmd;
		p_cent.Zfi = e_fi_max - (idx + 0.5)*lmd;
		Axis_slx p_cent_axis = dRadiustoAxis_gpu_slx(p_cent, direction_radius);
		squares[id].CenterRay.x = p_cent_axis.p[0];
		squares[id].CenterRay.y = p_cent_axis.p[1];
		squares[id].CenterRay.z = p_cent_axis.p[2];
	}

	id = idy * (width + 1) + idx;
	if (idx < width + 1 && idy < height + 1)//声线管束的四个角顶射线
	{
		p_corner.Yst = e_st_min + idy * lmd;
		p_corner.Zfi = e_fi_max - idx * lmd;
		Axis_slx p_corner_axis = dRadiustoAxis_gpu_slx(p_corner, direction_radius);
		rays[id].p.x = p_corner_axis.p[0];
		rays[id].p.y = p_corner_axis.p[1];
		rays[id].p.z = p_corner_axis.p[2];
	}
}

void create_virtualface_gpu(Direction *rays2, Square *squares2, int dwidth, int dheight, float lmd, Radius direction_radius,
	float e_st_min, float e_fi_max, int device_id)
{
	dim3 gridsize((dwidth + 1) / TILE_DIM + 1, (dheight + 1) / TILE_DIM + 1);
	dim3 blocksize(TILE_DIM, TILE_DIM);
	dcreate_virtualface_gpu << <gridsize, blocksize >> >(direction_radius, rays2, squares2, dheight, dwidth, e_st_min, e_fi_max, lmd);
	HANDLE_ERROR(cudaGetLastError());
}

__global__ void dStFiExtreme(Radius* d_spherical, float* d_st_min, float* d_st_max, float* d_fi_min, float* d_fi_max, int nodeNum)
{
	__shared__ float st_min_data[512];
	__shared__ float st_max_data[512];
	__shared__ float fi_min_data[512];
	__shared__ float fi_max_data[512];

	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	int tid = threadIdx.x;
	if (idx < nodeNum)
	{
		st_min_data[tid] = d_spherical[idx].Yst;
		st_max_data[tid] = d_spherical[idx].Yst;
		fi_min_data[tid] = d_spherical[idx].Zfi;
		fi_max_data[tid] = d_spherical[idx].Zfi;
	}
	__syncthreads();

	for (int s = blockDim.x / 2; s > 0; s >>= 1)
	{
		if (tid < s){
			st_min_data[tid] = (st_min_data[tid] < st_min_data[tid + s]) ? st_min_data[tid] : st_min_data[tid + s];
			st_max_data[tid] = (st_max_data[tid] > st_max_data[tid + s]) ? st_max_data[tid] : st_max_data[tid + s];
			fi_min_data[tid] = (fi_min_data[tid] < fi_min_data[tid + s]) ? fi_min_data[tid] : fi_min_data[tid + s];
			fi_max_data[tid] = (fi_max_data[tid] > fi_max_data[tid + s]) ? fi_max_data[tid] : fi_max_data[tid + s];
		}
		__syncthreads();
	}

	if (tid == 0)
	{
		d_st_min[blockIdx.x] = st_min_data[0];
		d_st_max[blockIdx.x] = st_max_data[0];
		d_fi_min[blockIdx.x] = fi_min_data[0];
		d_fi_max[blockIdx.x] = fi_max_data[0];
	}
}

__global__ void dStFiExtreme2(float* d_st_min, float* d_st_max, float* d_fi_min, float* d_fi_max, int nodeNum)
{
	__shared__ float st_min_data[512];
	__shared__ float st_max_data[512];
	__shared__ float fi_min_data[512];
	__shared__ float fi_max_data[512];

	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	int tid = threadIdx.x;
	if (idx < nodeNum)
	{
		st_min_data[tid] = d_st_min[idx];
		st_max_data[tid] = d_st_max[idx];
		fi_min_data[tid] = d_fi_min[idx];
		fi_max_data[tid] = d_fi_max[idx];
	}
	__syncthreads();

	for (int s = blockDim.x / 2; s > 0; s >>= 1)
	{
		if (tid < s){
			st_min_data[tid] = (st_min_data[tid] < st_min_data[tid + s]) ? st_min_data[tid] : st_min_data[tid + s];
			st_max_data[tid] = (st_max_data[tid] > st_max_data[tid + s]) ? st_max_data[tid] : st_max_data[tid + s];
			fi_min_data[tid] = (fi_min_data[tid] < fi_min_data[tid + s]) ? fi_min_data[tid] : fi_min_data[tid + s];
			fi_max_data[tid] = (fi_max_data[tid] > fi_max_data[tid + s]) ? fi_max_data[tid] : fi_max_data[tid + s];
		}
		__syncthreads();
	}

	if (tid == 0)
	{
		d_st_min[blockIdx.x] = st_min_data[0];
		d_st_max[blockIdx.x] = st_max_data[0];
		d_fi_min[blockIdx.x] = fi_min_data[0];
		d_fi_max[blockIdx.x] = fi_max_data[0];
	}
}
void getWidthHeight(float far_dis,Element *points, int num, int *width_calc, int *height_calc,int st, float fi, float lmd_calc, int *dheight_calc, float *de_st_min, float *de_fi_max, int device_num, int* divided_num, int max_pipeline_capicity)
{
	Radius direction = { far_dis, st, fi };	//射线表示用角度（r，st，fi）  θ =st，和z轴的夹角；为从正z轴来看自x轴按逆时针方向转到OM所转过的角  // snw 不能接立即数
	
	Radius* d_direction;
	HANDLE_ERROR(cudaMalloc((void**)&d_direction, sizeof(Radius)));
	HANDLE_ERROR(cudaMemcpy(d_direction, &direction, sizeof(Radius), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaDeviceSynchronize());

	dim3 threadSize(512, 1, 1);
	dim3 blockSize(num / 512 + 1, 1, 1);

	int numItem = 0;
	Element* d_point;
	HANDLE_ERROR(cudaMalloc((void**)&d_point, num * sizeof(Element)));
	//HANDLE_ERROR(cudaMemset(d_point, 0, num * sizeof(Element)));
	HANDLE_ERROR(cudaMemcpy(d_point, points, num * sizeof(Element), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaDeviceSynchronize());

	Radius* d_spherical;//读取数据时，球面坐标系矢量
	HANDLE_ERROR(cudaMalloc((void**)&d_spherical, num * sizeof(Radius)));
	HANDLE_ERROR(cudaMemset(d_spherical, 0, num * sizeof(Radius)));

	//将模型各点的坐标转换到球坐标系
	dAxistoRadius_gpu << <blockSize, threadSize >> >(d_spherical, d_point, direction, num);
	HANDLE_ERROR(cudaGetLastError());
	HANDLE_ERROR(cudaDeviceSynchronize());

	float* d_st_min;
	float* d_fi_max;
	float* d_st_max;
	float* d_fi_min;
	HANDLE_ERROR(cudaMalloc((void**)&d_st_min, blockSize.x * sizeof(float)));
	HANDLE_ERROR(cudaMalloc((void**)&d_fi_max, blockSize.x * sizeof(float)));
	HANDLE_ERROR(cudaMalloc((void**)&d_st_max, blockSize.x * sizeof(float)));
	HANDLE_ERROR(cudaMalloc((void**)&d_fi_min, blockSize.x * sizeof(float)));
	HANDLE_ERROR(cudaMemset(d_st_min, 0, blockSize.x * sizeof(float)));
	HANDLE_ERROR(cudaMemset(d_fi_max, 0, blockSize.x * sizeof(float)));
	HANDLE_ERROR(cudaMemset(d_st_max, 0, blockSize.x * sizeof(float)));
	HANDLE_ERROR(cudaMemset(d_fi_min, 0, blockSize.x * sizeof(float)));

	//两级归约求st,fi的最大最小值
	dStFiExtreme << <blockSize, threadSize >> >(d_spherical, d_st_min, d_st_max, d_fi_min, d_fi_max, num);
	HANDLE_ERROR(cudaGetLastError());
	HANDLE_ERROR(cudaDeviceSynchronize());
	numItem = blockSize.x;
	blockSize.x = numItem / threadSize.x + 1;
	while (numItem > 1)
	{
		dStFiExtreme2 << <blockSize, threadSize >> >(d_st_min, d_st_max, d_fi_min, d_fi_max, numItem);
		cudaDeviceSynchronize();
		numItem = blockSize.x;
		blockSize.x = numItem / threadSize.x + 1;
	}

	//读取三角面元的位置信息，转化到虚拟孔径面找到最值
	float e_st_min, e_st_max, e_fi_min, e_fi_max;
	HANDLE_ERROR(cudaMemcpy(&e_st_min, d_st_min, sizeof(float), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(&e_st_max, d_st_max, sizeof(float), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(&e_fi_min, d_fi_min, sizeof(float), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(&e_fi_max, d_fi_max, sizeof(float), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaDeviceSynchronize());
	*width_calc = (e_fi_max - e_fi_min + lmd_calc - 0.001) / lmd_calc;
	*height_calc = (e_st_max - e_st_min + lmd_calc - 0.001) / lmd_calc;
	*dheight_calc = (*height_calc + device_num - 0.001) / device_num;
	//printf("e_fi_max:%f, e_fi_min:%f",e_fi_max, e_fi_min);
	//printf("e_st_max:%f, e_st_min:%f",e_st_max, e_st_min);
	*de_st_min = e_st_min;
	*de_fi_max = e_fi_max;

	*divided_num = ceil((float)*width_calc * (float)*height_calc / (float)max_pipeline_capicity);
	

	cudaFree(d_direction);
	cudaFree(d_point);
	cudaFree(d_spherical);
	cudaFree(d_st_min);
	cudaFree(d_fi_max);
	cudaFree(d_st_max);
	cudaFree(d_fi_min);
}

void getWidthHeight_Sun(float far_dis,Box b, int *width, int *height, int st, float fi, float lmd, int *dheight, float *de_st_min, float *de_fi_max, int device_num)
{
	float  distance = far_dis;//snw
	float p[3];
	float dir[3];
	p[2] = distance*cos(st*PI/180);
	p[1] = distance*sin(st*PI / 180)*cos(fi*PI / 180);
	p[0] = distance*sin(st*PI / 180)*sin(fi*PI / 180);

	float s = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
	dir[0] = -p[0] / s;
	dir[1] = -p[1] / s;
	dir[2] = -p[2] / s;
	float dirx[3], diry[3];
	dirx[0] = dir[1] / (sqrt(dir[0] * dir[0] + dir[1] * dir[1]));
	dirx[1] = -dir[0] / (sqrt(dir[0] * dir[0] + dir[1] * dir[1]));
	dirx[2] = 0;

	diry[0] = dirx[1] * dir[2] - dirx[2] * dir[1];
	diry[1] = dirx[2] * dir[0] - dirx[0] * dir[2];
	diry[2] = dirx[0] * dir[1] - dirx[1] * dir[0];


	float D = -(dir[0] * p[0]) - (dir[1] * p[1]) - (dir[2] * p[2]);

	float point[8][3];
	point[0][0] = b.bmin[0]; point[0][1] = b.bmin[1], point[0][2] = b.bmin[2];
	point[1][0] = b.bmax[0]; point[1][1] = b.bmin[1], point[1][2] = b.bmin[2];
	point[2][0] = b.bmin[0]; point[2][1] = b.bmax[1], point[2][2] = b.bmin[2];
	point[3][0] = b.bmin[0]; point[3][1] = b.bmin[1], point[3][2] = b.bmax[2];
	point[4][0] = b.bmax[0]; point[4][1] = b.bmax[1], point[4][2] = b.bmin[2];
	point[5][0] = b.bmax[0]; point[5][1] = b.bmin[1], point[5][2] = b.bmax[2];
	point[6][0] = b.bmin[0]; point[6][1] = b.bmax[1], point[6][2] = b.bmax[2];
	point[7][0] = b.bmax[0]; point[7][1] = b.bmax[1], point[7][2] = b.bmax[2];

	float xmin=999999999, xmax=-999999999, ymin=9999999999, ymax=-9999999999;
	for (int i = 0; i < 8; i++)
	{
		float t = (dir[0] * point[i][0] + dir[1] * point[i][1] + dir[2] * point[i][2] + D)
			/ (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
		float dir_t[3];
		float x = point[i][0] - t*dir[0];
		float y = point[i][1] - t*dir[1];
		float z = point[i][2] - t*dir[2];

		dir_t[0] = x - p[0];
		dir_t[1] = y - p[1];
		dir_t[2] = z - p[2];

		float yy = dir_t[2] / diry[2];
		float xx = (dir_t[0] - y*diry[0]) / dirx[0];
		if (xx > xmax)
			xmax = xx;
		if (xx < xmin)
			xmin = xx;
		if (yy > ymax)
			ymax = yy;
		if (yy < ymin)
			ymin = yy;
	}
	*height = (ymax - ymin) / lmd;
	*width = (xmax - xmin) / lmd;
	*dheight = (*height + device_num - 1) / device_num;
	*de_st_min = ymin;
	*de_fi_max = xmax;

	return;
}

void ConstructCore(DynamicPlane *node, float e_st_min, float e_fi_max, int width, int height)
{
	node->st_min = e_st_min;
	node->fi_max = e_fi_max;
	node->height = height;
	node->width = width;
}

void ConstructNodeLeft(DynamicPlane parent, DynamicPlane* child, BinaryTimeTree* pre_angle_time, float lmd)
{
	if (parent.width >= parent.height)
	{
		child->flag = 0; // 列向划分

		int width;
		if (pre_angle_time->runtime == 0 || pre_angle_time->flag == 1)
		{
			// 第一度仿真或与前一度划分方向不一样时
			child->GPUnum = ceil(parent.GPUnum / 2); // ceil:返回大于等于它的最小整数
			child->ratio = ceil(parent.GPUnum / 2) / parent.GPUnum;
			width = ceil(parent.width * child->ratio);
		}
		else
		{
			int deltaW = (int)(0.5 * (pre_angle_time->rightchild->runtime * pre_angle_time->leftchild->GPUnum - pre_angle_time->leftchild->runtime * pre_angle_time->rightchild->GPUnum) / \
				(pre_angle_time->AveTimeLeft * pre_angle_time->rightchild->GPUnum + pre_angle_time->AveTimeRight * pre_angle_time->leftchild->GPUnum));
			width = (int)(parent.width * pre_angle_time->ratio) + deltaW;
			child->ratio = (float)width / parent.width;
		}
		float fi_max = parent.fi_max - (parent.width - width) * lmd;
		ConstructCore(child, parent.st_min, fi_max, width, parent.height);
	}
	else
	{
		child->flag = 1; // 行向划分

		int height;
		if (pre_angle_time->runtime == 0 || pre_angle_time->flag == 0)
		{
			// 第一度仿真或与前一度划分方向不一样时
			child->GPUnum = ceil(parent.GPUnum / 2); // ceil:返回大于等于它的最小整数
			child->ratio = ceil(parent.GPUnum / 2) / parent.GPUnum;
			height = ceil(parent.height * child->ratio);
		}
		else
		{
			int deltaH = (int)(0.5 * (pre_angle_time->rightchild->runtime * pre_angle_time->leftchild->GPUnum - pre_angle_time->leftchild->runtime * pre_angle_time->rightchild->GPUnum) / \
				(pre_angle_time->AveTimeLeft * pre_angle_time->rightchild->GPUnum + pre_angle_time->AveTimeRight * pre_angle_time->leftchild->GPUnum));
			height = (int)(parent.height*pre_angle_time->ratio) + deltaH;
			child->ratio = (float)height / parent.height;
		}
		float st_min = parent.st_min + (parent.height - height) * lmd;
		ConstructCore(child, st_min, parent.fi_max, parent.width, height);
	}
}

void ConstructNodeRight(DynamicPlane parent, DynamicPlane* child, BinaryTimeTree* pre_angle_time)
{
	if (parent.width >= parent.height)
	{
		child->flag = 0; // 按列

		int width;
		if (pre_angle_time->runtime == 0 || pre_angle_time->flag == 1)
		{
			// 第一度仿真或与前一度划分方向不一样时
			child->GPUnum = floor(parent.GPUnum / 2); // floor:返回小于等于它的最小整数
			width = parent.width - ceil(parent.width * ceil(parent.GPUnum / 2) / parent.GPUnum);
		}
		else
		{
			int deltaW = (int)(0.5 * (pre_angle_time->leftchild->runtime * pre_angle_time->rightchild->GPUnum - pre_angle_time->rightchild->runtime * pre_angle_time->leftchild->GPUnum) / \
				(pre_angle_time->AveTimeLeft * pre_angle_time->rightchild->GPUnum + pre_angle_time->AveTimeRight * pre_angle_time->leftchild->GPUnum));
			width = parent.width - (int)(parent.width * pre_angle_time->ratio) + deltaW;
		}
		ConstructCore(child, parent.st_min, parent.fi_max, width, parent.height);
	}
	else
	{
		child->flag = 1; // 按行

		int height;
		if (pre_angle_time->runtime == 0 || pre_angle_time->flag == 0)
		{
			child->GPUnum = floor(parent.GPUnum / 2); // floor:返回小于等于它的最小整数
			height = parent.height - ceil(parent.height * ceil(parent.GPUnum / 2) / parent.GPUnum);
		}
		else
		{
			int deltaH = (int)(0.5 * (pre_angle_time->leftchild->runtime * pre_angle_time->rightchild->GPUnum - pre_angle_time->rightchild->runtime * pre_angle_time->leftchild->GPUnum) / \
				(pre_angle_time->AveTimeLeft * pre_angle_time->rightchild->GPUnum + pre_angle_time->AveTimeRight * pre_angle_time->leftchild->GPUnum));
			height = parent.height - (int)(parent.height * pre_angle_time->ratio) + deltaH;
		}
		ConstructCore(child, parent.st_min, parent.fi_max, parent.width, height);
	}
}

/**************************
名称：DynamicPlane* ConstructVirtualFace()
描述：动态生成子孔径面
参数：BinaryTimeTree** pre_angle_time:前一角度的计算时间; int DeviceCount:GPU数量; float e_st_min, float e_fi_max, int width, int height:虚拟孔径面信息; float lmd:划分步长
返回值：DynamicPlane* dData:各GPU卡上的子孔径面边界信息
***************************/
void ConstructVirtualFace(DynamicPlane* array, DynamicPlane* dData, BinaryTimeTree** pre_angle_time, int DeviceCount, float e_st_min, float e_fi_max, int width, int height, float lmd)
{
	int NodeNum = 2 * DeviceCount - 1;

	//DynamicPlane* array = (DynamicPlane*)malloc(NodeNum * sizeof(DynamicPlane));
	int index = 1;

	while (index <= NodeNum)
	{
		if (index == 1) // 创建根节点
		{
			array[index - 1].st_min = e_st_min;
			array[index - 1].fi_max = e_fi_max;
			array[index - 1].height = height;
			array[index - 1].width = width;
			array[index - 1].GPUnum = DeviceCount;
			
			if (index >= DeviceCount && index <= NodeNum) // 单卡
			{
				dData[index - DeviceCount] = array[index - 1];
			}
			index++;
			continue;
		}

		if (index % 2 == 0) // 左子节点
		{
			ConstructNodeLeft(array[index / 2 - 1], &array[index - 1], pre_angle_time[index / 2 - 1], lmd);
		}
		else // 右子节点
		{
			ConstructNodeRight(array[index / 2 - 1], &array[index - 1], pre_angle_time[index / 2 - 1]);
		}

		if (index >= DeviceCount && index <= NodeNum)
		{
			dData[index - DeviceCount] = array[index - 1];
		}
		index++;
	}
}

/**************************
名称：BinaryTimeTree** ConstructTimeTree()
描述：将各卡的计算时间存储成二叉树结构
参数：DynamicPlane* plane:子孔径面信息; float *runtime:各卡计算时间; int DeviceCount:GPU数量
返回值：BinaryTimeTree** timetree:前一角度的各卡计算时间
***************************/
BinaryTimeTree** ConstructTimeTree(DynamicPlane* plane, float *runtime, int DeviceCount)
{
	int NodeNum = 2 * DeviceCount - 1;
	BinaryTimeTree** timetree = (BinaryTimeTree**)malloc(NodeNum * sizeof(BinaryTimeTree*));
	int i = 1;
	while (i <= NodeNum)
	{
		if (i == 1)
		{
			BinaryTimeTree *root = (BinaryTimeTree*)malloc(sizeof(BinaryTimeTree));
			root->runtime = 0;
			timetree[i - 1] = root;
			i++;
			continue;
		}
		timetree[i - 1] = (BinaryTimeTree*)malloc(sizeof(BinaryTimeTree));
		BinaryTimeTree *temp = new BinaryTimeTree();
		if (i % 2 == 0)
			(timetree[i / 2 - 1])->leftchild = temp;
		else
			(timetree[i / 2 - 1])->rightchild = temp;
		timetree[i - 1] = temp;

		if (i >= DeviceCount && i <= NodeNum)
		{
			(timetree[i - 1])->runtime = runtime[i - DeviceCount];
			(timetree[i - 1])->GPUnum = 1;
		}
		i++;
	}
	i = NodeNum - 1;
	while (i > 0)
	{
		(timetree[i / 2 - 1])->flag = plane[i].flag;
		(timetree[i / 2 - 1])->runtime = (timetree[i])->runtime + (timetree[i - 1])->runtime;
		(timetree[i / 2 - 1])->ratio = plane[i - 1].ratio;
		(timetree[i / 2 - 1])->GPUnum = (timetree[i])->GPUnum + (timetree[i - 1])->GPUnum;
		(timetree[i / 2 - 1])->AveTimeRight = (timetree[i])->runtime / (plane[i].flag == 0 ? plane[i].width : plane[i].height);
		(timetree[i / 2 - 1])->AveTimeLeft = (timetree[i - 1])->runtime / (plane[i - 1].flag == 0 ? plane[i - 1].width : plane[i - 1].height);
		i -= 2;
	}
	return timetree;
}