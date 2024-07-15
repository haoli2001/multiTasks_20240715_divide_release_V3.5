#include<cuda_runtime.h>
#include<device_launch_parameters.h>
#include<device_functions.h>
#include<string.h>
#include "integral_gpu.h"
#include "handlerror.h"
#include "sm_20_atomic_functions.h"

__global__ void reduce_add_re(float* d_sum_re, ReimOutput* d_in, int raysBeamNum)
{
	__shared__ float sdata[512];

	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	int tid = threadIdx.x;

	sdata[tid] = (idx < raysBeamNum) ? d_in[idx].re : 0;
	__syncthreads();

	for (int s = blockDim.x / 2; s > 0; s >>= 1)
	{
		if (tid < s){
			sdata[tid] +=  sdata[tid + s];
		}
		__syncthreads();
	}

	if (tid == 0)
	{
		d_sum_re[blockIdx.x] = sdata[0];
	}
}

__global__ void reduce_re(float* d_sum_re, float* d_in, int Num)
{
	__shared__ float sdata[512];

	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	int tid = threadIdx.x;

	sdata[tid] = (idx < Num) ? d_in[idx] : 0;
	__syncthreads();

	for (int s = blockDim.x / 2; s > 0; s >>= 1)
	{
		if (tid < s){
			sdata[tid] += sdata[tid + s];
		}
		__syncthreads();
	}

	if (tid == 0)
	{
		d_sum_re[blockIdx.x] = sdata[0];
	}
}

__global__ void reduce_add_im(float* d_sum_im, ReimOutput* d_in, int raysBeamNum)
{
	__shared__ float sdata[512];

	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	int tid = threadIdx.x;

	sdata[tid] = (idx < raysBeamNum) ? d_in[idx].im : 0;
	__syncthreads();

	for (int s = blockDim.x / 2; s > 0; s >>= 1)
	{
		if (tid < s){
			sdata[tid] += sdata[tid + s];
		}
		__syncthreads();
	}

	if (tid == 0)
	{
		d_sum_im[blockIdx.x] = sdata[0];
	}
}

__global__ void reduce_im(float* d_sum_im, float* d_in, int Num)
{
	__shared__ float sdata[512];

	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	int tid = threadIdx.x;

	sdata[tid] = (idx < Num) ? d_in[idx] : 0;
	__syncthreads();

	for (int s = blockDim.x / 2; s > 0; s >>= 1)
	{
		if (tid < s){
			sdata[tid] += sdata[tid + s];
		}
		__syncthreads();
	}

	if (tid == 0)
	{
		d_sum_im[blockIdx.x] = sdata[0];
	}
}

//�������
comp product(comp c1, comp c2)
{
	comp x = {0,0};
	x.re = c1.re * c2.re - c1.im * c2.im;
	x.im = c1.im * c2.re + c2.im * c1.re;
	return x;
}

__device__ comp complex_product(comp c1, comp c2)
{
	comp x = {0,0};
	x.re = c1.re * c2.re - c1.im * c2.im;
	x.im = c1.im * c2.re + c2.im * c1.re;
	return x;
}

__device__ comp complex_add(comp c1, comp c2)
{
	comp sum = {0,0};
	sum.re = c1.re + c2.re;
	sum.im = c1.im + c2.im;
	return sum;
}

__device__ float sa_function(float f1, float f2, float wavenum)
{
	float result = 0;
	if (f2 == 0)
	{
		result = f1;
	}
	else
	{
		result = f1 * (sinf(wavenum / 2 * f2) / (wavenum / 2 * f2));
	}
	return result;
}

__global__ void integral(RayBeamInfo* d_effrays, int raysBeamNum, float wavenumber, ReimOutput* d_reim)
{
	float WaveNumber = 2 * wavenumber;
	float d_omgdelt1, d_omgdelt2, tmp;
	comp d_ctemp, d_reim_tmp = {0,0};

	PointCoor point0, point1, point2;
	Vector normal_vector, center_vecter;//向量求面积的参数，声束管线中心点位置
	PointCoor launch_point, nomal_launch_point;//声源点位置，归一化之后声源点位置20201014
	float launch_point_model;
	
	for (int idx = threadIdx.x + blockIdx.x*blockDim.x; idx < raysBeamNum; idx += blockDim.x*gridDim.x)
	{
		if(d_effrays[idx].p_cent_distance == 0)//20190116
			return;
			
		PointCoor3 parameter = d_effrays[idx].parameter;
		float p_cent_distance = d_effrays[idx].p_cent_distance;
		PointCoor p_2d[5];
		p_2d[0] = d_effrays[idx].point_2D[0];
		p_2d[1] = d_effrays[idx].point_2D[1];
		p_2d[2] = d_effrays[idx].point_2D[2];
		p_2d[3] = d_effrays[idx].point_2D[3];
		p_2d[4] = d_effrays[idx].point_2D[0];
			
		for (int i = 0; i < 4; i++)
		{
			d_omgdelt1 = parameter.p[0] * (p_2d[i + 1].p[0] - p_2d[i].p[0]) + parameter.p[1] * (p_2d[i + 1].p[1] - p_2d[i].p[1]);
			d_omgdelt2 = parameter.p[1] * (p_2d[i + 1].p[0] - p_2d[i].p[0]) - parameter.p[0] * (p_2d[i + 1].p[1] - p_2d[i].p[1]);
			tmp = sa_function(d_omgdelt2, d_omgdelt1, WaveNumber);
			d_ctemp.re = 0;
			d_ctemp.im = -(parameter.p[0] * WaveNumber * (p_2d[i + 1].p[0] + p_2d[i].p[0]) / 2 + parameter.p[1] * WaveNumber * (p_2d[i + 1].p[1] + p_2d[i].p[1]) / 2);
			d_ctemp.re = cosf(d_ctemp.im) * tmp;
			d_ctemp.im = sinf(d_ctemp.im) * tmp;
			d_reim_tmp = complex_add(d_reim_tmp, d_ctemp);
		}
		parameter.p[0] = (parameter.p[0] == 0) ? 9999999999 : parameter.p[0];
		parameter.p[1] = (parameter.p[1] == 0) ? 9999999999 : parameter.p[1];

		d_reim_tmp.re = d_reim_tmp.re / WaveNumber / (parameter.p[0] * parameter.p[0] + parameter.p[1] * parameter.p[1]);
		d_reim_tmp.im = d_reim_tmp.im / WaveNumber / (parameter.p[0] * parameter.p[0] + parameter.p[1] * parameter.p[1]);

		d_ctemp.re = 0;
		d_ctemp.im = 2 * wavenumber * p_cent_distance;
		d_ctemp.re = cosf(d_ctemp.im);
		d_ctemp.im = sinf(d_ctemp.im);

		d_reim_tmp = complex_product(d_reim_tmp, d_ctemp);
		d_ctemp.re = -1 / p_cent_distance / p_cent_distance / p_cent_distance / (-2 * PI);
		d_ctemp.im = wavenumber / p_cent_distance / p_cent_distance / (-2 * PI);
		d_reim_tmp = complex_product(d_reim_tmp, d_ctemp);

		d_reim[idx].re = d_reim_tmp.re * parameter.p[2];
		d_reim[idx].im = d_reim_tmp.im * parameter.p[2];
		//2023.2.9
		d_effrays[idx].gorden_re = d_reim[idx].re;
		d_effrays[idx].gorden_im = d_reim[idx].im;

		launch_point = d_effrays[idx].launch_point;
		launch_point_model = sqrt(launch_point.p[0] * launch_point.p[0] + launch_point.p[1] * launch_point.p[1] + launch_point.p[2] * launch_point.p[2]);
		nomal_launch_point.p[0] = launch_point.p[0] / launch_point_model;
		nomal_launch_point.p[1] = launch_point.p[1] / launch_point_model;
		nomal_launch_point.p[2] = launch_point.p[2] / launch_point_model;
		center_vecter.x = 1.0f / 4.0f * (d_effrays[idx].ray_index[0].p[0] + d_effrays[idx].ray_index[1].p[0] + d_effrays[idx].ray_index[2].p[0] + d_effrays[idx].ray_index[3].p[0]);
		center_vecter.y = 1.0f / 4.0f * (d_effrays[idx].ray_index[0].p[1] + d_effrays[idx].ray_index[1].p[1] + d_effrays[idx].ray_index[2].p[1] + d_effrays[idx].ray_index[3].p[1]);
		center_vecter.z = 1.0f / 4.0f * (d_effrays[idx].ray_index[0].p[2] + d_effrays[idx].ray_index[1].p[2] + d_effrays[idx].ray_index[2].p[2] + d_effrays[idx].ray_index[3].p[2]);
		d_effrays[idx].Z0 = center_vecter.x * nomal_launch_point.p[0] + center_vecter.y * nomal_launch_point.p[1] + center_vecter.z * nomal_launch_point.p[2];


		point0 = d_effrays[idx].ray_index[0];
		point1 = d_effrays[idx].ray_index[1];
		point2 = d_effrays[idx].ray_index[2];
		normal_vector.x = (point0.p[1] - point2.p[1]) * (point1.p[2] - point0.p[2]) - (point1.p[1] - point0.p[1]) * (point0.p[2] - point2.p[2]);
		normal_vector.y = (point0.p[2] - point2.p[2]) * (point1.p[0] - point0.p[0]) - (point1.p[2] - point0.p[2]) * (point0.p[0] - point2.p[0]);
		normal_vector.z = (point0.p[0] - point2.p[0]) * (point1.p[1] - point0.p[1]) - (point1.p[0] - point0.p[0]) * (point0.p[1] - point2.p[1]);
		d_effrays[idx].gorden = sqrt(normal_vector.x * normal_vector.x + normal_vector.y * normal_vector.y + normal_vector.z * normal_vector.z);
	}
}

__global__ void parameter(MatStruct* transMat, RayBeamInfo* d_effrays, int raysBeamNum)
{
	float tmp_data;
	
	for (int idx = threadIdx.x + blockIdx.x*blockDim.x; idx < raysBeamNum; idx += blockDim.x*gridDim.x)
	{
		if(d_effrays[idx].p_cent_distance == 0)
			return;
			
		MatStruct resg_mat = transMat[idx];
		PointCoor tmp_ray = d_effrays[idx].launch_point;
		float dis = d_effrays[idx].p_cent_distance;
		for (int k = 0; k < 3; k++)
		{
			tmp_data = 0;
			for (int i = 0; i < 4; i++)
			{
				tmp_data += resg_mat.p[(4 * k) + i] * tmp_ray.p[i];
			}
			d_effrays[idx].parameter.p[k] = tmp_data / dis;
		}
	}
}

__global__ void coordinate2D(MatStruct* transMat, RayBeamInfo* d_effrays, int raysBeamNum)
{
	float tmp_data;
	PointCoor tmp_ray;
	__shared__ MatStruct s_mat[8];

	int tid = threadIdx.x;
	int bdx = blockIdx.x; 

	if ((bdx * 8 + tid / 4) >= raysBeamNum || d_effrays[bdx * 8 + tid / 4].p_cent_distance == 0)
		return;

	s_mat[tid / 4] = transMat[bdx * 8 + tid / 4];
	__syncthreads();

	tmp_ray = d_effrays[bdx * 8 + tid / 4].ray_index[tid % 4];
	for (int k = 0; k < 3; k++)
	{
		tmp_data = 0;
		for (int i = 0; i < 4; i++)
		{
			tmp_data += s_mat[tid/4].p[(4 * k) + i] * tmp_ray.p[i];
		}
		d_effrays[bdx * 8 + tid / 4].point_2D[tid % 4].p[k] = tmp_data;
	}
	
	
	/*MatStruct resg_mat;//20190427
	float tmp_data;
	PointCoor tmp_ray;
	for (int idx = threadIdx.x + blockIdx.x*blockDim.x; idx < raysBeamNum; idx += blockDim.x*gridDim.x)
	{
		//if (d_effrays[idx].p_cent_distance == 0)//20190116
		//	return;
			
		resg_mat = transMat[idx];

		for (int j = 0; j < 4; j++)
		{
			tmp_ray = d_effrays[idx].ray_index[j];
			for (int k = 0; k < 3; k++)
			{
				tmp_data = 0;
				for (int i = 0; i < 4; i++)
				{
					tmp_data += resg_mat.p[(4 * k) + i] * tmp_ray.p[i];
				}
				d_effrays[idx].point_2D[j].p[k] = tmp_data;
			}
		}
	}*/

}


__global__ void build_transMat(RayBeamInfo* rays, int raysBeamNum, Vector* d_axis, Vector* d_center, MatStruct* transMat)
{

	//20191218

	float tmp;
	Vector vector_x, vector_y, vector_z;
	Vector normal_vector;

	PointCoor point0, point1, point2;
	Vector point_center;
	
	for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < raysBeamNum; idx += gridDim.x * blockDim.x)//2018-12-11
	{

		if(rays[idx].p_cent_distance == 0)//20190116
			return;

		point0 = rays[idx].ray_index[0];
		point1 = rays[idx].ray_index[1];
		point2 = rays[idx].ray_index[2];

		normal_vector.x = (point0.p[1] - point2.p[1]) * (point1.p[2] - point0.p[2]) - (point1.p[1] - point0.p[1]) * (point0.p[2] - point2.p[2]);
		normal_vector.y = (point0.p[2] - point2.p[2]) * (point1.p[0] - point0.p[0]) - (point1.p[2] - point0.p[2]) * (point0.p[0] - point2.p[0]);
		normal_vector.z = (point0.p[0] - point2.p[0]) * (point1.p[1] - point0.p[1]) - (point1.p[0] - point0.p[0]) * (point0.p[1] - point2.p[1]);
		tmp = normal_vector.x * normal_vector.x + normal_vector.y * normal_vector.y + normal_vector.z * normal_vector.z;
		vector_z.x = normal_vector.x / sqrtf(tmp);
		vector_z.y = normal_vector.y / sqrtf(tmp);
		vector_z.z = normal_vector.z / sqrtf(tmp);//z

		point_center = d_center[idx];

		tmp = (point0.p[0] - point_center.x)*(point0.p[0] - point_center.x) \
			+ (point0.p[1] - point_center.y)*(point0.p[1] - point_center.y) \
			+ (point0.p[2] - point_center.z)*(point0.p[2] - point_center.z);

		//ͶӰ����ϵx��
		vector_x.x = (point0.p[0] - point_center.x) / sqrtf(tmp);
		vector_x.y = (point0.p[1] - point_center.y) / sqrtf(tmp);
		vector_x.z = (point0.p[2] - point_center.z) / sqrtf(tmp);//x

		//ͶӰ����ϵy��
		vector_y.x = -(vector_x.y * vector_z.z - vector_z.y * vector_x.z);
		vector_y.y = -(vector_x.z * vector_z.x - vector_z.z * vector_x.x);
		vector_y.z = -(vector_x.x * vector_z.y - vector_z.x * vector_x.y);//y

		transMat[idx].p[0] = vector_x.x;//x
		transMat[idx].p[1] = vector_x.y;
		transMat[idx].p[2] = vector_x.z;
		transMat[idx].p[3] = -point_center.x * vector_x.x - point_center.y * vector_x.y - point_center.z * vector_x.z;
		transMat[idx].p[4] = vector_y.x;//y
		transMat[idx].p[5] = vector_y.y;
		transMat[idx].p[6] = vector_y.z;
		transMat[idx].p[7] = -point_center.x * vector_y.x - point_center.y * vector_y.y - point_center.z * vector_y.z;
		transMat[idx].p[8] = vector_z.x;//z
		transMat[idx].p[9] = vector_z.y;
		transMat[idx].p[10] = vector_z.z;
		transMat[idx].p[11] = -point_center.x * vector_z.x - point_center.y * vector_z.y - point_center.z * vector_z.z;
	}
}

__global__ void copy_data_gpu(Direction *d_rays, Square *d_squares, RayBeamInfo* d_rays_save, int raysBeamNum, Vector* d_center, PointCoor launchpoint_tmp, ReimOutput* d_reim, ConfigStruct* config)
{
	float CSpeed = 1500.0;//shengsu 

	//float fs = config->sampling_frequency;//caiyanglv 姬梓遇20210831
	for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < raysBeamNum; idx += gridDim.x * blockDim.x)//2018-12-11
	{
		//if ((d_squares[idx].right == true) && (d_rays[d_squares[idx].ray_index[4]].triangle_index != -1))
		if ((d_squares[idx].right == true))
		{
			float ray00 = d_rays[d_squares[idx].CornerRayIndex.x].p.x;
			float ray01 = d_rays[d_squares[idx].CornerRayIndex.x].p.y;
			float ray02 = d_rays[d_squares[idx].CornerRayIndex.x].p.z;
			int triangle_index = d_rays[d_squares[idx].CornerRayIndex.x].triangle_index;
			float ray10 = d_rays[d_squares[idx].CornerRayIndex.y].p.x;
			float ray11 = d_rays[d_squares[idx].CornerRayIndex.y].p.y;
			float ray12 = d_rays[d_squares[idx].CornerRayIndex.y].p.z;
			float ray20 = d_rays[d_squares[idx].CornerRayIndex.z].p.x;
			float ray21 = d_rays[d_squares[idx].CornerRayIndex.z].p.y;
			float ray22 = d_rays[d_squares[idx].CornerRayIndex.z].p.z;
			float ray30 = d_rays[d_squares[idx].CornerRayIndex.w].p.x;
			float ray31 = d_rays[d_squares[idx].CornerRayIndex.w].p.y;
			float ray32 = d_rays[d_squares[idx].CornerRayIndex.w].p.z;
			d_center[idx].x = 1.0f / 4.0f * (ray00 + ray10 + ray20 + ray30);
			d_center[idx].y = 1.0f / 4.0f * (ray01 + ray11 + ray21 + ray31);
			d_center[idx].z = 1.0f / 4.0f * (ray02 + ray12 + ray22 + ray32);
			d_rays_save[idx].p_cent_distance = sqrt(pow(d_center[idx].x - launchpoint_tmp.p[0], 2) + pow(d_center[idx].y - launchpoint_tmp.p[1], 2) + pow(d_center[idx].z - launchpoint_tmp.p[2], 2));
			d_rays_save[idx].launch_point = launchpoint_tmp;
			d_rays_save[idx].ray_index[0].p[0] = ray00;
			d_rays_save[idx].ray_index[0].p[1] = ray01;
			d_rays_save[idx].ray_index[0].p[2] = ray02;
			d_rays_save[idx].ray_index[0].p[3] = 1.0f;
			d_rays_save[idx].ray_index[1].p[0] = ray10;
			d_rays_save[idx].ray_index[1].p[1] = ray11;
			d_rays_save[idx].ray_index[1].p[2] = ray12;
			d_rays_save[idx].ray_index[1].p[3] = 1.0f;
			d_rays_save[idx].ray_index[2].p[0] = ray20;
			d_rays_save[idx].ray_index[2].p[1] = ray21;
			d_rays_save[idx].ray_index[2].p[2] = ray22;
			d_rays_save[idx].ray_index[2].p[3] = 1.0f;
			d_rays_save[idx].ray_index[3].p[0] = ray30;
			d_rays_save[idx].ray_index[3].p[1] = ray31;
			d_rays_save[idx].ray_index[3].p[2] = ray32;
			d_rays_save[idx].ray_index[3].p[3] = 1.0f;
			d_reim[idx].triangle_index = triangle_index;
			//d_rays_save[idx].launch_point.p[0] = x;//20190423
			//d_rays_save[idx].launch_point.p[1] = y;
			//d_rays_save[idx].launch_point.p[2] = z;
			//d_rays_save[idx].launch_point.p[3] = 1.0f;
			/*d_rays_save[idx].launch_point = launchpoint_tmp;
			d_rays_save[idx].ray_index[0].p[0] = d_rays[d_squares[idx].ray_index[0]].p[0];
			d_rays_save[idx].ray_index[0].p[1] = d_rays[d_squares[idx].ray_index[0]].p[1];
			d_rays_save[idx].ray_index[0].p[2] = d_rays[d_squares[idx].ray_index[0]].p[2];
			d_rays_save[idx].ray_index[0].p[3] = 1.0f;
			d_rays_save[idx].ray_index[1].p[0] = d_rays[d_squares[idx].ray_index[1]].p[0];
			d_rays_save[idx].ray_index[1].p[1] = d_rays[d_squares[idx].ray_index[1]].p[1];
			d_rays_save[idx].ray_index[1].p[2] = d_rays[d_squares[idx].ray_index[1]].p[2];
			d_rays_save[idx].ray_index[1].p[3] = 1.0f;
			d_rays_save[idx].ray_index[2].p[0] = d_rays[d_squares[idx].ray_index[2]].p[0];
			d_rays_save[idx].ray_index[2].p[1] = d_rays[d_squares[idx].ray_index[2]].p[1];
			d_rays_save[idx].ray_index[2].p[2] = d_rays[d_squares[idx].ray_index[2]].p[2];
			d_rays_save[idx].ray_index[2].p[3] = 1.0f;
			d_rays_save[idx].ray_index[3].p[0] = d_rays[d_squares[idx].ray_index[3]].p[0];
			d_rays_save[idx].ray_index[3].p[1] = d_rays[d_squares[idx].ray_index[3]].p[1];
			d_rays_save[idx].ray_index[3].p[2] = d_rays[d_squares[idx].ray_index[3]].p[2];
			d_rays_save[idx].ray_index[3].p[3] = 1.0f;
			d_center[idx].x = 1.0f / 4.0f * (d_rays_save[idx].ray_index[0].p[0] + d_rays_save[idx].ray_index[1].p[0] + d_rays_save[idx].ray_index[2].p[0] + d_rays_save[idx].ray_index[3].p[0]);
			d_center[idx].y = 1.0f / 4.0f * (d_rays_save[idx].ray_index[0].p[1] + d_rays_save[idx].ray_index[1].p[1] + d_rays_save[idx].ray_index[2].p[1] + d_rays_save[idx].ray_index[3].p[1]);
			d_center[idx].z = 1.0f / 4.0f * (d_rays_save[idx].ray_index[0].p[2] + d_rays_save[idx].ray_index[1].p[2] + d_rays_save[idx].ray_index[2].p[2] + d_rays_save[idx].ray_index[3].p[2]);
			d_rays_save[idx].p_cent_distance = sqrt(pow(d_center[idx].x - launchpoint_tmp.p[0], 2) + pow(d_center[idx].y - launchpoint_tmp.p[1], 2) + pow(d_center[idx].z - launchpoint_tmp.p[2], 2));*/
		}
	}
}

comp sound_field_integral_gpu(Direction *d_rays, Square *d_squares, float wavelength, int raysBeamNum, RayBeamInfo* d_effrays, Vector* d_center, Vector* d_axis, 
	MatStruct* d_transMat, ReimOutput* d_reim, float* d_sum_re, float* d_sum_im, float fi, float si,  ConfigStruct* config)
{
	float wavenumber = (2 * PI) / wavelength;

	
	float x = config->far_distance * sin(PI / 180 * fi)*cos(PI / 180 * si);
	float y = config->far_distance * sin(PI / 180 * fi)*sin(PI / 180 * si);
	float z = config->far_distance * cos(PI / 180 * fi);
	PointCoor launchpoint_tmp = {x, y, z, 1.0f};

	dim3 threadSize(16, 1, 1);
	dim3 blockSize(raysBeamNum / 16 + 1, 1, 1);
	
	//��ʱ
	//cudaEvent_t start, stop;
	//HANDLE_ERROR(cudaEventCreate(&start));
	//HANDLE_ERROR(cudaEventCreate(&stop));
	//HANDLE_ERROR(cudaEventRecord(start, 0));
	//���ݿ�������

	copy_data_gpu << <blockSize, threadSize >> >(d_rays, d_squares, d_effrays, raysBeamNum, d_center, launchpoint_tmp, d_reim, config);


	//HANDLE_ERROR(cudaDeviceSynchronize());
	//HANDLE_ERROR(cudaEventRecord(stop, 0));
	//HANDLE_ERROR(cudaEventSynchronize(stop));
	//float elapsedTime;
	//HANDLE_ERROR(cudaEventElapsedTime(&elapsedTime, start, stop));
	//printf("Time to copy data at %d: %fms\n", degree, elapsedTime);
	//HANDLE_ERROR(cudaEventDestroy(start));
	//HANDLE_ERROR(cudaEventDestroy(stop));

	build_transMat << <blockSize, threadSize >> >(d_effrays, raysBeamNum, d_axis, d_center, d_transMat);//20191218
	//HANDLE_ERROR(cudaDeviceSynchronize());

	//compute_axis_x_gpu << <blockSize, threadSize >> >(d_effrays, raysBeamNum, d_center, d_axis);
	//HANDLE_ERROR(cudaDeviceSynchronize());

	//build_transMat << <blockSize, threadSize >> >(d_effrays, d_transMat, raysBeamNum, d_center, d_axis);
	//HANDLE_ERROR(cudaDeviceSynchronize());

	dim3 threadSize1(32, 1, 1);
	dim3 blockSize1(raysBeamNum/8+1, 1, 1);
	coordinate2D << <blockSize1, threadSize1 >> >(d_transMat, d_effrays, raysBeamNum);
	//HANDLE_ERROR(cudaDeviceSynchronize());
	
	//coordinate2D << <raysBeamNum / 64 + 1, 64 >> >(d_transMat, d_effrays, raysBeamNum);

	parameter << <blockSize, threadSize >> >(d_transMat, d_effrays, raysBeamNum);
	//HANDLE_ERROR(cudaDeviceSynchronize());

	integral << <blockSize, threadSize >> >(d_effrays, raysBeamNum, wavenumber, d_reim);
	//HANDLE_ERROR(cudaDeviceSynchronize());

	int numItem = 0;
	
	dim3 THREADSIZE(512, 1, 1);
	dim3 BLOCKSIZE(raysBeamNum / 512 + 1, 1, 1);
	//ReduceAdd.
	reduce_add_re<<<BLOCKSIZE,THREADSIZE>>>(d_sum_re, d_reim, raysBeamNum);
	cudaDeviceSynchronize();
	reduce_add_im<<<BLOCKSIZE,THREADSIZE>>>(d_sum_im, d_reim, raysBeamNum);
	cudaDeviceSynchronize();
	numItem = BLOCKSIZE.x;
	BLOCKSIZE.x = numItem / THREADSIZE.x + 1;

	while (numItem > 1)
	{
		reduce_re<<<BLOCKSIZE,THREADSIZE>>>(d_sum_re, d_sum_re, numItem);
		cudaDeviceSynchronize();
		reduce_im<<<BLOCKSIZE,THREADSIZE>>>(d_sum_im, d_sum_im, numItem);
		cudaDeviceSynchronize();
		numItem = BLOCKSIZE.x;
		BLOCKSIZE.x = numItem / THREADSIZE.x + 1;
	}

	float sum_re = 0;
	float sum_im = 0;

     

	HANDLE_ERROR(cudaMemcpy(&sum_re, d_sum_re, sizeof(float), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(&sum_im, d_sum_im, sizeof(float), cudaMemcpyDeviceToHost));



	//����ǿ�ȼ���
	comp sum = { sum_re, sum_im };

	return sum;
}

comp sound_field_integral_gpu_DivRay(Direction *d_rays, Square *d_squares, float wavelength, int raysBeamNum,
	RayBeamInfo* d_effrays, Vector* d_center, Vector* d_axis, MatStruct* d_transMat, ReimOutput* d_reim, float* d_sum_re, float* d_sum_im, float fi, float si,  ConfigStruct* config)
{
	raysBeamNum *= 4;
	
	float x = config->far_distance * sin(PI / 180 * fi)*cos(PI / 180 * si);
	float y = config->far_distance * sin(PI / 180 * fi)*sin(PI / 180 * si);
	float z = config->far_distance * cos(PI / 180 * fi);
	PointCoor launchpoint_tmp = {x, y, z, 1.0f};
	
	float wavenumber = (2 * PI) / wavelength;

	dim3 threadSize(16, 1, 1);
	dim3 blockSize(raysBeamNum / 16 + 1, 1, 1);

	copy_data_gpu<< <blockSize, threadSize >> >(d_rays, d_squares, d_effrays, raysBeamNum, d_center, launchpoint_tmp, d_reim, config);
	//HANDLE_ERROR(cudaDeviceSynchronize());

	/*RayBeamInfo * h_effrays;
	h_effrays = (RayBeamInfo*)malloc(15976009 * sizeof(RayBeamInfo));
	HANDLE_ERROR(cudaMemcpy(h_effrays, d_effrays, 15976009 * sizeof(RayBeamInfo), cudaMemcpyDeviceToHost));
	
	printf("�ĸ��Ƕ�����---------------------\n");
	for(int m=1050; m<1100;m++)
	{
		printf("%d ",m);
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				printf("%.7f ", h_effrays[m].ray_index[i].p[j]);
			}
			printf("\n");
		}
	}
	printf("---------------------\n");*/

	build_transMat << <blockSize, threadSize >> >(d_effrays, raysBeamNum, d_axis, d_center, d_transMat);
	//HANDLE_ERROR(cudaDeviceSynchronize());

	//compute_axis_x_gpu << <blockSize, threadSize >> >(d_effrays, raysBeamNum, d_center, d_axis);
	//HANDLE_ERROR(cudaDeviceSynchronize());
	
	//build_transMat << <blockSize, threadSize >> >(d_effrays, d_transMat, raysBeamNum, d_center, d_axis);
	//HANDLE_ERROR(cudaDeviceSynchronize());

	//float* h_transMat;
	//h_transMat = (float*)malloc(320356 * 12 * sizeof(float));
	//HANDLE_ERROR(cudaMemcpy2D(h_transMat, 12 * sizeof(float), d_transMat, pitch, 12 * sizeof(float), 320356, cudaMemcpyDeviceToHost));
	//
	//printf("�任����---------------------\n");
	//for (int i = 164; i < 165; i++)
	//{
	//	for (int j = 0; j < 12; ++j)
	//	{
	//		printf("%.7f ", h_transMat[i*(12) + j]);
	//		if ((j == 3) || (j == 7) || (j == 11))
	//			printf("\n");
	//	}
	//}
	//printf("---------------------\n");

	//dim3 ThreadPerBlock(1, 1, 1);
	//dim3 BlockPerGrid(raysBeamNum, 1, 1);
	//coordinate2D << <blockSize, threadSize >> >(d_transMat, d_effrays, raysBeamNum);
	dim3 threadSize1(32, 1, 1);
	dim3 blockSize1(raysBeamNum / 8 + 1, 1, 1);
	coordinate2D << <blockSize1, threadSize1 >> >(d_transMat, d_effrays, raysBeamNum);
	//HANDLE_ERROR(cudaDeviceSynchronize());

	//coordinate2D << <blockSize, threadSize >> >(d_transMat, d_effrays, raysBeamNum);

	parameter << <blockSize, threadSize >> >(d_transMat, d_effrays, raysBeamNum);
	//HANDLE_ERROR(cudaDeviceSynchronize());

	/*HANDLE_ERROR(cudaMemcpy(h_effrays, d_effrays, 15976009 * sizeof(RayBeamInfo), cudaMemcpyDeviceToHost));
	printf("u,v,w = ");
	for(int m=1050; m<1100;m++)
	{
		printf("%d ",m);
		for (int j = 0; j < 3; j++)
		{
			printf("%.7f ", h_effrays[m].parameter.p[j]);
		}
		printf("\n");
	}
	printf("\n---------------------\n");*/

	integral << <blockSize, threadSize >> >(d_effrays, raysBeamNum, wavenumber, d_reim);
	//HANDLE_ERROR(cudaDeviceSynchronize());
	
	/*float sum_reim_re = 0;
	float sum_reim_im = 0;
	comp* h_reim;
	h_reim = (comp*)malloc(15976009 * sizeof(comp));
	HANDLE_ERROR(cudaMemcpy(h_reim, d_reim, 15976009 * sizeof(comp), cudaMemcpyDeviceToHost));
	for (int i = 0; i < raysBeamNum; i++)
	{
		sum_reim_re += h_reim[i].re;
		sum_reim_im += h_reim[i].im;
		printf("%d reim = %.16f,%.16f\n",i, h_reim[i].re, h_reim[i].im);
	}
	printf("reim = %.16f,%.16f\n", sum_reim_re, sum_reim_im) ;*/
	
	int numItem = 0;
	dim3 THREADSIZE(512, 1, 1);
	dim3 BLOCKSIZE(raysBeamNum / 512 + 1, 1, 1);
	//ReduceAdd.
	reduce_add_re << <BLOCKSIZE, THREADSIZE >> >(d_sum_re, d_reim, raysBeamNum);
	cudaDeviceSynchronize();
	reduce_add_im << <BLOCKSIZE, THREADSIZE >> >(d_sum_im, d_reim, raysBeamNum);
	cudaDeviceSynchronize();
	numItem = BLOCKSIZE.x;
	BLOCKSIZE.x = numItem / THREADSIZE.x + 1;

	while (numItem > 1)
	{
		reduce_re << <BLOCKSIZE, THREADSIZE >> >(d_sum_re, d_sum_re, numItem);
		cudaDeviceSynchronize();
		reduce_im << <BLOCKSIZE, THREADSIZE >> >(d_sum_im, d_sum_im, numItem);
		cudaDeviceSynchronize();
		numItem = BLOCKSIZE.x;
		BLOCKSIZE.x = numItem / THREADSIZE.x + 1;
	}

	float sum_re = 0;
	float sum_im = 0;
	HANDLE_ERROR(cudaMemcpy(&sum_re, d_sum_re, sizeof(float), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(&sum_im, d_sum_im, sizeof(float), cudaMemcpyDeviceToHost));

	//����ǿ�ȼ���
	comp sum = { sum_re, sum_im };
	//printf("integral1: %.16f,%.16f\n", sum.re, sum.im);
	return sum;
}

comp sound_field_integral_gpu_DivRay3(Direction *d_rays, Square *d_squares, float wavelength, int raysBeamNum,
	RayBeamInfo* d_effrays, Vector* d_center, Vector* d_axis, MatStruct* d_transMat, ReimOutput* d_reim, float* d_sum_re, float* d_sum_im, float fi, float si,  ConfigStruct* config)
{
	raysBeamNum *= 4;
	
	float x = config->far_distance * sin(PI / 180 * fi)*cos(PI / 180 * si);  //snw
	float y = config->far_distance * sin(PI / 180 * fi)*sin(PI / 180 * si);
	float z = config->far_distance * cos(PI / 180 * fi);
	PointCoor launchpoint_tmp = {x, y, z, 1.0f};
	
	float wavenumber = (2 * PI) / wavelength;

	dim3 threadSize(16, 1, 1);
	dim3 blockSize(raysBeamNum / 16 + 1, 1, 1);

	copy_data_gpu<< <blockSize, threadSize >> >(d_rays, d_squares, d_effrays, raysBeamNum, d_center, launchpoint_tmp, d_reim, config);
	//HANDLE_ERROR(cudaDeviceSynchronize());

	/*RayBeamInfo * h_effrays;
	h_effrays = (RayBeamInfo*)malloc(15976009 * sizeof(RayBeamInfo));
	HANDLE_ERROR(cudaMemcpy(h_effrays, d_effrays, 15976009 * sizeof(RayBeamInfo), cudaMemcpyDeviceToHost));
	
	printf("�ĸ��Ƕ�����---------------------\n");
	for(int m=9716; m<9717;m++)
	{
		//printf("%d ",m);
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				printf("%.7f ", h_effrays[m].ray_index[i].p[j]);
			}
			printf("\n");
		}
	}
	printf("---------------------\n");
	printf("distance:%f\n",h_effrays[9716].p_cent_distance);*/
	
	/*Vector * h_center;
	h_center = (Vector*)malloc(15976009 * sizeof(Vector));
	HANDLE_ERROR(cudaMemcpy(h_center, d_center, 15976009 * sizeof(Vector), cudaMemcpyDeviceToHost));
	printf("���ĵ�����---------------------\n");
	for(int m=9716; m<9717;m++)
	{
		printf("%.7f %.7f %.7f\n", h_center[m].x, h_center[m].y, h_center[m].z);
	}
	printf("---------------------\n");*/

	build_transMat << <blockSize, threadSize >> >(d_effrays, raysBeamNum, d_axis, d_center, d_transMat);
	//HANDLE_ERROR(cudaDeviceSynchronize());
	
	/*Vector * h_axis;
	h_axis = (Vector*)malloc(3 * 15976009 * sizeof(Vector));
	HANDLE_ERROR(cudaMemcpy(h_axis, d_axis, 3 * 15976009 * sizeof(Vector), cudaMemcpyDeviceToHost));
	printf("---------------------\n");
	for(int m=9716; m<9717;m++)
	{
		printf("%.7f %.7f %.7f\n", h_axis[m].x, h_axis[m].y, h_axis[m].z);
	}
	printf("---------------------\n");*/

	//compute_axis_x_gpu << <blockSize, threadSize >> >(d_effrays, raysBeamNum, d_center, d_axis);
	//HANDLE_ERROR(cudaDeviceSynchronize());

	/*HANDLE_ERROR(cudaMemcpy(h_axis, d_axis, 3 * 15976009 * sizeof(Vector), cudaMemcpyDeviceToHost));
	printf("---------------------\n");
	for(int m=9716+9564644*2; m<9717+9564644*2;m++)
	{
		printf("%.7f %.7f %.7f\n", h_axis[m].x, h_axis[m].y, h_axis[m].z);
	}
	for(int m=9716+9564644; m<9717+9564644;m++)
	{
		printf("%.7f %.7f %.7f\n", h_axis[m].x, h_axis[m].y, h_axis[m].z);
	}
	printf("---------------------\n");*/
	
	//build_transMat << <blockSize, threadSize >> >(d_effrays, d_transMat, raysBeamNum, d_center, d_axis);
	//HANDLE_ERROR(cudaDeviceSynchronize());

	/*MatStruct* h_transMat;
	h_transMat = (MatStruct*)malloc(15976009 * sizeof(MatStruct));
	HANDLE_ERROR(cudaMemcpy(h_transMat, d_transMat, sizeof(MatStruct) * 15976009, cudaMemcpyDeviceToHost));
	
	printf("�任����---------------------\n");
	for (int i = 9716; i < 9717; i++)
	{
		for (int j = 0; j < 12; ++j)
		{
			printf("%.7f ", h_transMat[i].p[j]);
			if ((j == 3) || (j == 7) || (j == 11))
				printf("\n");
		}
	}
	printf("---------------------\n");*/

	//dim3 ThreadPerBlock(1, 1, 1);
	//dim3 BlockPerGrid(raysBeamNum, 1, 1);
	//coordinate2D << <blockSize, threadSize >> >(d_transMat, d_effrays, raysBeamNum);
	dim3 threadSize1(32, 1, 1);
	dim3 blockSize1(raysBeamNum / 8 + 1, 1, 1);
	coordinate2D << <blockSize1, threadSize1 >> >(d_transMat, d_effrays, raysBeamNum);
	//HANDLE_ERROR(cudaDeviceSynchronize());
	
	//coordinate2D << <blockSize, threadSize >> >(d_transMat, d_effrays, raysBeamNum);

	parameter << <blockSize, threadSize >> >(d_transMat, d_effrays, raysBeamNum);
	//HANDLE_ERROR(cudaDeviceSynchronize());

	/*HANDLE_ERROR(cudaMemcpy(h_effrays, d_effrays, 15976009 * sizeof(RayBeamInfo), cudaMemcpyDeviceToHost));
	printf("u,v,w = ");
	for(int m=9716; m<9717;m++)
	{
		printf("%d ",m);
		for (int j = 0; j < 3; j++)
		{
			printf("%.7f ", h_effrays[m].parameter.p[j]);
		}
		printf("\n");
	}
	printf("\n---------------------\n");*/

	integral << <blockSize, threadSize >> >(d_effrays, raysBeamNum, wavenumber, d_reim);
	//HANDLE_ERROR(cudaDeviceSynchronize());
	
	/*float sum_reim_re = 0;
	float sum_reim_im = 0;
	comp* h_reim;
	h_reim = (comp*)malloc(15976009 * sizeof(comp));
	HANDLE_ERROR(cudaMemcpy(h_reim, d_reim, 15976009 * sizeof(comp), cudaMemcpyDeviceToHost));
	for (int i = 9716; i < 9717; i++)
	{
		//sum_reim_re += h_reim[i].re;
		//sum_reim_im += h_reim[i].im;
		printf("%d reim = %.16f,%.16f\n",i, h_reim[i].re, h_reim[i].im);
	}
	//printf("reim = %.16f,%.16f\n", sum_reim_re, sum_reim_im) ;*/
	
	int numItem = 0;
	dim3 THREADSIZE(512, 1, 1);
	dim3 BLOCKSIZE(raysBeamNum / 512 + 1, 1, 1);
	//ReduceAdd.
	reduce_add_re << <BLOCKSIZE, THREADSIZE >> >(d_sum_re, d_reim, raysBeamNum);
	cudaDeviceSynchronize();
	reduce_add_im << <BLOCKSIZE, THREADSIZE >> >(d_sum_im, d_reim, raysBeamNum);
	cudaDeviceSynchronize();
	numItem = BLOCKSIZE.x;
	BLOCKSIZE.x = numItem / THREADSIZE.x + 1;

	while (numItem > 1)
	{
		reduce_re << <BLOCKSIZE, THREADSIZE >> >(d_sum_re, d_sum_re, numItem);
		cudaDeviceSynchronize();
		reduce_im << <BLOCKSIZE, THREADSIZE >> >(d_sum_im, d_sum_im, numItem);
		cudaDeviceSynchronize();
		numItem = BLOCKSIZE.x;
		BLOCKSIZE.x = numItem / THREADSIZE.x + 1;
	}

	float sum_re = 0;
	float sum_im = 0;
	HANDLE_ERROR(cudaMemcpy(&sum_re, d_sum_re, sizeof(float), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(&sum_im, d_sum_im, sizeof(float), cudaMemcpyDeviceToHost));

	//����ǿ�ȼ���
	comp sum = { sum_re, sum_im };
	//printf("integral1: %.16f,%.16f\n", sum.re, sum.im);
	return sum;
}

float TS_compute(float far_dis,comp sum, float wavelength)  //����ǿ�ȼ���
{
	comp reim = { 0, 0 };
	float mag = 0;
	float intensity = 0;
	reim.re = 0;
	reim.im = -(2 * PI) / wavelength * far_dis;
	reim.re = cosf(reim.im);
	reim.im = sinf(reim.im);
	sum = product(sum, reim);
	mag = sqrt(sum.re * sum.re + sum.im * sum.im);
	if (mag < 99999999999.0)
		intensity = 20 * log10f(mag * far_dis * far_dis);
	return intensity;
}

//float sound_field_integral_gpu_f(Direction *d_rays, Square *d_squares, float wavelength, int width, int height)
//{
//	//Direction   * d_rays;
//	//Square      * d_squares;
//	RayBeamInfo * d_effrays;
//	Vector      * d_center;
//	Vector      * d_axis_z;
//	PointCoor   * ac;
//	PointCoor   * bd;
//	float       * nv_length;
//	Vector      * normal_vector;
//	Vector      * d_axis_x;
//	float       * d_vectorX_length;
//	Vector      * d_axis_y;
//	size_t      pitch;
//	size_t      pitch1;
//	size_t 	    pitch2;
//	float       * d_transMat;
//	float intensity = 0;
//	int nrow = 4;
//	int ncol = 4;
//	float wavenumber = (2 * PI) / wavelength;
//	int num = 0;
//	int raysNum = (width + 1) * (height + 1) + width * height;
//	int raysBeamNum = width * height;
//
//	/*HANDLE_ERROR(cudaMalloc((void**)&d_rays, raysNum * sizeof(Direction)));
//	HANDLE_ERROR(cudaMalloc((void**)&d_squares, raysBeamNum * sizeof(Square)));*/
//	HANDLE_ERROR(cudaMalloc((void**)&d_effrays, raysBeamNum * sizeof(RayBeamInfo)));
//	HANDLE_ERROR(cudaMalloc((void**)&d_center, raysBeamNum * sizeof(Vector)));
//	HANDLE_ERROR(cudaMalloc((void**)&d_axis_z, raysBeamNum * sizeof(Vector)));
//	HANDLE_ERROR(cudaMalloc((void**)&ac, raysBeamNum * sizeof(PointCoor)));
//	HANDLE_ERROR(cudaMalloc((void**)&bd, raysBeamNum * sizeof(PointCoor)));
//	HANDLE_ERROR(cudaMalloc((void**)&nv_length, raysBeamNum * sizeof(float)));
//	HANDLE_ERROR(cudaMalloc((void**)&normal_vector, raysBeamNum * sizeof(Vector)));
//	HANDLE_ERROR(cudaMalloc((void**)&d_axis_x, raysBeamNum * sizeof(Vector)));
//	HANDLE_ERROR(cudaMalloc((void**)&d_vectorX_length, raysBeamNum * sizeof(float)));
//	HANDLE_ERROR(cudaMalloc((void**)&d_axis_y, raysBeamNum * sizeof(Vector)));
//	HANDLE_ERROR(cudaMallocPitch((void**)&d_transMat, &pitch, nrow * ncol * sizeof(float), raysBeamNum));
//
//	/*HANDLE_ERROR(cudaMemcpy(d_rays, rays, raysNum * sizeof(Direction), cudaMemcpyHostToDevice));
//	HANDLE_ERROR(cudaMemcpy(d_squares, squares, raysBeamNum * sizeof(Square), cudaMemcpyHostToDevice));*/
//
//	HANDLE_ERROR(cudaMemset(d_effrays, 0, raysBeamNum * sizeof(RayBeamInfo)));
//	HANDLE_ERROR(cudaMemset(d_center, 0, raysBeamNum * sizeof(Vector)));
//	HANDLE_ERROR(cudaMemset(d_axis_z, 0, raysBeamNum * sizeof(Vector)));
//	HANDLE_ERROR(cudaMemset(ac, 0, raysBeamNum * sizeof(PointCoor)));
//	HANDLE_ERROR(cudaMemset(bd, 0, raysBeamNum * sizeof(PointCoor)));
//	HANDLE_ERROR(cudaMemset(nv_length, 0, raysBeamNum * sizeof(float)));
//	HANDLE_ERROR(cudaMemset(normal_vector, 0, raysBeamNum * sizeof(Vector)));
//	HANDLE_ERROR(cudaMemset(d_axis_x, 0, raysBeamNum * sizeof(Vector)));
//	HANDLE_ERROR(cudaMemset(d_vectorX_length, 0, raysBeamNum * sizeof(float)));
//	HANDLE_ERROR(cudaMemset(d_axis_y, 0, raysBeamNum * sizeof(Vector)));
//	HANDLE_ERROR(cudaMemset2D(d_transMat, pitch, 0, nrow * ncol * sizeof(float), raysBeamNum));
//
//	dim3 threadSize(512, 1, 1);
//	dim3 blockSize(width * height / 512 + 1, 1, 1);
//
//	//��ʱ
//	//cudaEvent_t start, stop;
//	//HANDLE_ERROR(cudaEventCreate(&start));
//	//HANDLE_ERROR(cudaEventCreate(&stop));
//	//HANDLE_ERROR(cudaEventRecord(start, 0));
//	//���ݿ�������
//	copy_data_gpu << <blockSize, threadSize >> >(d_rays, d_squares, d_effrays, width, height);
//	cudaDeviceSynchronize();
//
//
//	//HANDLE_ERROR(cudaEventRecord(stop, 0));
//	//HANDLE_ERROR(cudaEventSynchronize(stop));
//	//float elapsedTime;
//	//HANDLE_ERROR(cudaEventElapsedTime(&elapsedTime, start, stop));
//	//printf("Time to copy data at %d: %fms\n", degree, elapsedTime);
//	//HANDLE_ERROR(cudaEventDestroy(start));
//	//HANDLE_ERROR(cudaEventDestroy(stop));
//
//	//RayBeamInfo * h_effrays;
//	//h_effrays = (RayBeamInfo*)malloc(raysBeamNum*sizeof(RayBeamInfo));
//	//HANDLE_ERROR(cudaMemcpy(h_effrays, d_effrays, raysBeamNum*sizeof(RayBeamInfo), cudaMemcpyDeviceToHost));
//
//	//printf("�ĸ��Ƕ�����---------------------\n");
//	//for (int i = 0; i < 4; i++)
//	//{
//	//for (int j = 0; j < 4; j++)
//	//{
//	//printf("%.7f ", h_effrays[30078].ray_index[i].p[j]);
//	//}
//	//printf("\n");
//	//}
//	//printf("---------------------\n");
//
//	compute_axis_z_gpu << <blockSize, threadSize >> >(d_effrays, width, height, ac, bd, nv_length, normal_vector, d_axis_z);
//	cudaDeviceSynchronize();
//	cudaFree(ac);
//	cudaFree(bd);
//	cudaFree(nv_length);
//	cudaFree(normal_vector);
//
//	compute_axis_x_gpu << <blockSize, threadSize >> >(d_effrays, width, height, d_center, d_vectorX_length, d_axis_x);
//	cudaDeviceSynchronize();
//
//	//Vector * h_center;
//	//h_center = (Vector*)malloc(raysBeamNum*sizeof(Vector));
//	//HANDLE_ERROR(cudaMemcpy(h_center, d_center, raysBeamNum*sizeof(Vector), cudaMemcpyDeviceToHost));
//	//printf("center = %.7f, %.7f, %.7f\n", h_center[30078].x, h_center[30078].y, h_center[30078].z);
//	//free(h_center);
//
//	//float * h_vectorX_length;
//	//h_vectorX_length = (float*)malloc(raysBeamNum*sizeof(float));
//	//HANDLE_ERROR(cudaMemcpy(h_vectorX_length, d_vectorX_length, raysBeamNum*sizeof(float), cudaMemcpyDeviceToHost));
//	//printf("length = %.7f\n", h_vectorX_length[30078]);
//	//free(h_vectorX_length);
//
//	cudaFree(d_vectorX_length);
//
//	compute_axis_y_gpu << <blockSize, threadSize >> >(d_axis_x, d_axis_z, width, height, d_axis_y);
//	cudaDeviceSynchronize();
//
//	build_transMat << <blockSize, threadSize >> >(d_transMat, width, height, d_axis_x, d_axis_y, d_axis_z, d_center, pitch);
//	cudaDeviceSynchronize();
//	//float* h_transMat;
//	//h_transMat = (float*)malloc(raysBeamNum * nrow * ncol * sizeof(float));
//	//HANDLE_ERROR(cudaMemcpy2D(h_transMat, nrow * ncol * sizeof(float), d_transMat, pitch, nrow * ncol * sizeof(float), raysBeamNum, cudaMemcpyDeviceToHost));
//
//	//printf("�任����---------------------\n");
//	//for (int i = 30078; i < 30079; i++)
//	//{
//	//for (int j = 0; j < 16; ++j)
//	//{
//	//printf("%.7f ", h_transMat[i*(nrow * ncol) + j]);
//	//if ((j == 3) || (j == 7) || (j == 11) || (j == 15))
//	//printf("\n");
//	//}
//	//}
//	//printf("---------------------\n");
//
//	cudaFree(d_axis_z);
//	cudaFree(d_axis_x);
//	cudaFree(d_axis_y);
//	cudaFree(d_center);
//
//	float* tempvar;
//	HANDLE_ERROR(cudaMallocPitch((void**)&tempvar, &pitch1, nrow * ncol * sizeof(float), raysBeamNum));
//	HANDLE_ERROR(cudaMemset2D(tempvar, pitch1, 0, nrow * ncol * sizeof(float), raysBeamNum));
//	coordinate2D << <blockSize, threadSize >> >(d_transMat, d_effrays, tempvar, raysBeamNum, pitch, pitch1);
//	cudaDeviceSynchronize();
//	cudaFree(tempvar);
//
//	distance << <blockSize, threadSize >> >(d_effrays, raysBeamNum);
//	cudaDeviceSynchronize();
//
//	float* tempVar;
//	HANDLE_ERROR(cudaMallocPitch((void**)&tempVar, &pitch2, ncol * sizeof(float), raysBeamNum));
//	HANDLE_ERROR(cudaMemset2D(tempVar, pitch2, 0, ncol * sizeof(float), raysBeamNum));
//	parameter << <blockSize, threadSize >> >(d_transMat, d_effrays, tempVar, raysBeamNum, pitch, pitch2);
//	cudaDeviceSynchronize();
//	cudaFree(d_transMat);
//	cudaFree(tempVar);
//
//	//HANDLE_ERROR(cudaMemcpy(h_effrays, d_effrays, raysBeamNum*sizeof(RayBeamInfo), cudaMemcpyDeviceToHost));
//	//printf("u,v,w = ");
//	//for (int j = 0; j < 3; j++)
//	//{
//	//printf("%.7f ", h_effrays[30078].parameter.p[j]);
//	//}
//	//printf("\n---------------------\n");
//
//	float       * d_omgdelt1;
//	float       * d_omgdelt2;
//	float       * d_tempVar;
//	comp        * d_ctemp;
//	comp        * d_reim;
//	HANDLE_ERROR(cudaMalloc((void**)&d_omgdelt1, raysBeamNum * sizeof(float)));
//	HANDLE_ERROR(cudaMalloc((void**)&d_omgdelt2, raysBeamNum * sizeof(float)));
//	HANDLE_ERROR(cudaMalloc((void**)&d_tempVar, raysBeamNum * sizeof(float)));
//	HANDLE_ERROR(cudaMalloc((void**)&d_ctemp, raysBeamNum * sizeof(comp)));
//	HANDLE_ERROR(cudaMalloc((void**)&d_reim, raysBeamNum * sizeof(comp)));
//	HANDLE_ERROR(cudaMemset(d_omgdelt1, 0, raysBeamNum * sizeof(float)));
//	HANDLE_ERROR(cudaMemset(d_omgdelt2, 0, raysBeamNum * sizeof(float)));
//	HANDLE_ERROR(cudaMemset(d_tempVar, 0, raysBeamNum * sizeof(float)));
//	HANDLE_ERROR(cudaMemset(d_ctemp, 0, raysBeamNum * sizeof(comp)));
//	HANDLE_ERROR(cudaMemset(d_reim, 0, raysBeamNum * sizeof(comp)));
//
//	integral << <blockSize, threadSize >> >(d_effrays, raysBeamNum, wavenumber, d_omgdelt1, d_omgdelt2, d_tempVar, d_ctemp, d_reim);
//	cudaDeviceSynchronize();
//
//	//comp* h_reim;
//	//h_reim = (comp*)malloc(raysBeamNum*sizeof(comp));
//	//HANDLE_ERROR(cudaMemcpy(h_reim, d_reim, raysBeamNum*sizeof(comp), cudaMemcpyDeviceToHost));
//	//for(int i=15000;i<15500;i++)
//	//{
//	//printf("%d reim = %.10f,%.10f\n",i, h_reim[i].re, h_reim[i].im);
//	//}
//
//	cudaFree(d_omgdelt1);
//	cudaFree(d_omgdelt2);
//	cudaFree(d_tempVar);
//	cudaFree(d_ctemp);
//
//	int numItem = 0;
//	float* d_sum_re;
//	float* d_sum_im;
//	HANDLE_ERROR(cudaMalloc((void**)&d_sum_re, blockSize.x * sizeof(float)));
//	HANDLE_ERROR(cudaMalloc((void**)&d_sum_im, blockSize.x * sizeof(float)));
//	HANDLE_ERROR(cudaMemset(d_sum_re, 0, blockSize.x * sizeof(float)));
//	HANDLE_ERROR(cudaMemset(d_sum_im, 0, blockSize.x * sizeof(float)));
//
//	//ReduceAdd.
//	reduce_add_re << <blockSize, threadSize >> >(d_sum_re, d_reim, raysBeamNum);
//	cudaDeviceSynchronize();
//	reduce_add_im << <blockSize, threadSize >> >(d_sum_im, d_reim, raysBeamNum);
//	cudaDeviceSynchronize();
//	numItem = blockSize.x;
//	blockSize.x = numItem / threadSize.x + 1;
//
//	while (numItem > 1)
//	{
//		reduce_re << <blockSize, threadSize >> >(d_sum_re, d_sum_re, numItem);
//		cudaDeviceSynchronize();
//		reduce_im << <blockSize, threadSize >> >(d_sum_im, d_sum_im, numItem);
//		cudaDeviceSynchronize();
//		numItem = blockSize.x;
//		blockSize.x = numItem / threadSize.x + 1;
//	}
//
//	float sum_re = 0;
//	float sum_im = 0;
//	HANDLE_ERROR(cudaMemcpy(&sum_re, d_sum_re, sizeof(float), cudaMemcpyDeviceToHost));
//	HANDLE_ERROR(cudaMemcpy(&sum_im, d_sum_im, sizeof(float), cudaMemcpyDeviceToHost));
//	//printf("sum_re = %.10f, sum_im = %.10f\n", sum_re, sum_im);
//
//	//����ǿ�ȼ���
//	comp Reim;
//	comp sum = { sum_re, sum_im };
//	float mag = 0;
//
//	Reim.re = 0;
//	Reim.im = -wavenumber * 1000;
//	Reim.re = cosf(Reim.im);
//	Reim.im = sinf(Reim.im);
//	//printf("Reim.re=%f, Reim.im=%f\n",Reim.re, Reim.im);
//	sum = product(sum, Reim);
//	mag = sqrt(sum.re * sum.re + sum.im * sum.im);
//	//printf("mag = %f\n",mag);
//	if (mag < 99999999999.0)
//		intensity = 20 * log10f(mag * 1000 * 1000);
//	//printf("The sound field intensity is %f\n", intensity);
//
//	//cudaFree(d_rays);
//	//cudaFree(d_squares);
//	cudaFree(d_effrays);
//	cudaFree(d_reim);
//	cudaFree(d_sum_re);
//	cudaFree(d_sum_im);
//	//free(h_transMat);
//	//free(h_effrays);
//	//free(h_reim);
//
//	return intensity;
//}
