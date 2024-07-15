#include<cuda_runtime.h>
#include<device_launch_parameters.h>
#include<device_functions.h>
#include "common_struct.h"
#include "virtualface_gpu.h"
#include "handlerror.h"
#include<stdio.h>


/*__global__ void CopyRayTube(Square *d_squares, int RayTubeNum, Square *d_squares2, Direction* d_rays2, Direction* d_rays, int* d_squares_pred, int* d_DivRayTubeNum, int RayTubeNumNew)
{
	int id = blockIdx.x * blockDim.x + threadIdx.x;

	//原声线管束中的射线，不需要重新追踪
	if (id < RayTubeNum && d_squares[id].IsDivRayTube == true)
	{
		int index = 4 * (d_DivRayTubeNum[blockIdx.x] + d_squares_pred[id]);

		int4 corner_ray_idx = d_squares[id].CornerRayIndex;
		int center_ray_idx = d_squares[id].CenterRayIndex;
		
		Direction ray0 = d_rays[corner_ray_idx.x];
		Direction ray1 = d_rays[corner_ray_idx.y];
		Direction ray2 = d_rays[corner_ray_idx.z];
		Direction ray3 = d_rays[corner_ray_idx.w];
		Direction ray4 = d_rays[center_ray_idx];
		
		d_rays2[index + RayTubeNumNew] = ray0;
		d_rays2[index + 1 + RayTubeNumNew] = ray1;
		d_rays2[index + 2 + RayTubeNumNew] = ray2;
		d_rays2[index + 3 + RayTubeNumNew] = ray3;
		d_rays2[index/4 + RayTubeNumNew*2] = ray4;
	}
}*/

__global__ void DivRayTubeCenter(Square *d_squares2, Direction* d_rays2, int RayTubeNumNew, int DivRayTubeNum, Axis direction)
{
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	int idx = id + RayTubeNumNew * 2 + DivRayTubeNum;

	//每根声线管束中心射线的坐标
	if (id < RayTubeNumNew)
	{
		//d_rays2[idx].caled = false;//20190427
		//d_rays2[idx].flag = false;
		//d_rays2[idx].distance = 0;
		//d_rays2[idx].times = 0;
		//d_rays2[idx].triangle_index = -1;
		//d_rays2[idx].dir[0] = -direction.x;
		//d_rays2[idx].dir[1] = -direction.y;
		//d_rays2[idx].dir[2] = -direction.z;

		int4 cr = d_squares2[id].CornerRayIndex;

		if ((id % 2) == 0)
		{
			float3 ray1 = d_rays2[cr.y].p;
			float3 ray3 = d_rays2[cr.w].p;

			d_squares2[id].CenterRay.x = 0.5 * (ray1.x + ray3.x);
			d_squares2[id].CenterRay.y = 0.5 * (ray1.y + ray3.y);
			d_squares2[id].CenterRay.z = 0.5 * (ray1.z + ray3.z);
		}
		if ((id % 2) == 1)
		{
			float3 ray0 = d_rays2[cr.x].p;
			float3 ray2 = d_rays2[cr.z].p;

			d_squares2[id].CenterRay.x = 0.5 * (ray0.x + ray2.x);
			d_squares2[id].CenterRay.y = 0.5 * (ray0.y + ray2.y);
			d_squares2[id].CenterRay.z = 0.5 * (ray0.z + ray2.z);
		}
	}
}

__global__ void CreateRayTubeInfo(Square *d_squares, Direction* d_rays, int RayTubeNum, Direction* d_rays2, int* d_squares_pred, 
	int* d_DivRayTubeNum, float lmd, float fi, Axis direction, int RayTubeNumNew)
{

	int id = blockIdx.x * blockDim.x + threadIdx.x;

	double ang = (double)fi;

	//新生成的需要重新追踪的射线
	if (id < RayTubeNum && d_squares[id].IsDivRayTube == true)
	{
		int index = 4 * (d_DivRayTubeNum[blockIdx.x] + d_squares_pred[id]);

		float3 CenterRay = d_squares[id].CenterRay;
		int4 corner_ray_idx = d_squares[id].CornerRayIndex;
		int center_ray_idx = d_squares[id].CenterRayIndex;

		Direction ray0 = d_rays[corner_ray_idx.x];
		Direction ray1 = d_rays[corner_ray_idx.y];
		Direction ray2 = d_rays[corner_ray_idx.z];
		Direction ray3 = d_rays[corner_ray_idx.w];
		Direction ray4 = d_rays[center_ray_idx];

		//新射线生成
		d_rays2[index].p.x = CenterRay.x;
		d_rays2[index].p.y = CenterRay.y;
		d_rays2[index].p.z = CenterRay.z + lmd;

		d_rays2[index + 1].p.x = CenterRay.x - lmd * sin(ang * 3.1415926 / 180);
		d_rays2[index + 1].p.y = CenterRay.y + lmd * cos(ang * 3.1415926 / 180);
		d_rays2[index + 1].p.z = CenterRay.z;

		d_rays2[index + 2].p.x = CenterRay.x;
		d_rays2[index + 2].p.y = CenterRay.y;
		d_rays2[index + 2].p.z = CenterRay.z - lmd;

		d_rays2[index + 3].p.x = CenterRay.x + lmd * sin(ang * 3.1415926 / 180);
		d_rays2[index + 3].p.y = CenterRay.y - lmd * cos(ang * 3.1415926 / 180);
		d_rays2[index + 3].p.z = CenterRay.z;

		//旧射线信息拷贝
		d_rays2[index + RayTubeNumNew] = ray0; // 20191223 将CopyRayTube和CreateRayTubeInfo合并为一个核函数
		d_rays2[index + 1 + RayTubeNumNew] = ray1;
		d_rays2[index + 2 + RayTubeNumNew] = ray2;
		d_rays2[index + 3 + RayTubeNumNew] = ray3;
		d_rays2[index / 4 + RayTubeNumNew * 2] = ray4;
	}
}

__global__ void CreateRayTubeIndex(Square *d_squares2, int RayTubeNumNew, int DivRayTubeNum)
{
	int id = blockIdx.x * blockDim.x + threadIdx.x;

	//每根声线管束的四个角顶射线编号
	if (id < RayTubeNumNew)
	{
		if ((id % 4) == 0)
		{
			d_squares2[id].CornerRayIndex.x = id + RayTubeNumNew;
			d_squares2[id].CornerRayIndex.y = id;
			d_squares2[id].CornerRayIndex.z = RayTubeNumNew * 2 + id / 4;
			d_squares2[id].CornerRayIndex.w = id + 3;
			d_squares2[id].CenterRayIndex = RayTubeNumNew * 2 + DivRayTubeNum + id;
		}
		if ((id % 4) == 1)
		{
			d_squares2[id].CornerRayIndex.x = id - 1;
			d_squares2[id].CornerRayIndex.y = id + RayTubeNumNew;
			d_squares2[id].CornerRayIndex.z = id;
			d_squares2[id].CornerRayIndex.w = RayTubeNumNew * 2 + id / 4;
			d_squares2[id].CenterRayIndex = RayTubeNumNew * 2 + DivRayTubeNum + id;
		}
		if ((id % 4) == 2)
		{
			d_squares2[id].CornerRayIndex.x = RayTubeNumNew * 2 + id / 4;
			d_squares2[id].CornerRayIndex.y = id - 1;
			d_squares2[id].CornerRayIndex.z = id + RayTubeNumNew;
			d_squares2[id].CornerRayIndex.w = id;
			d_squares2[id].CenterRayIndex = RayTubeNumNew * 2 + DivRayTubeNum + id;
		}
		if ((id % 4) == 3)
		{
			d_squares2[id].CornerRayIndex.x = id;
			d_squares2[id].CornerRayIndex.y = RayTubeNumNew * 2 + id / 4;
			d_squares2[id].CornerRayIndex.z = id - 1;
			d_squares2[id].CornerRayIndex.w = id + RayTubeNumNew;
			d_squares2[id].CenterRayIndex = RayTubeNumNew * 2 + DivRayTubeNum + id;
		}
	}
}

void DivRayTube(Direction *d_rays, Square *d_squares, Direction *d_rays2, Square *d_squares2, int* d_DivRayTubeNum, int DivRayTubeNum, int* d_sum_gmem, 
	int* d_squares_pred, int RayTubeNum, float lmd, float fi, Axis direction)
{
	int RayTubeNumNew = DivRayTubeNum * 4;
	//printf("RayTubeNumNew %d\n", RayTubeNumNew);
	
	dim3 SizeOfthread(512, 1, 1);
	dim3 SizeOfblock(RayTubeNumNew / 512 + 1, 1, 1);
	CreateRayTubeIndex << <SizeOfblock, SizeOfthread >> >(d_squares2, RayTubeNumNew, DivRayTubeNum);
	
	/*Square* squares2;
	squares2 = (Square*)malloc(3998000 * sizeof(Square));
	HANDLE_ERROR(cudaMemcpy(squares2, d_squares2, 3998000 * sizeof(Square), cudaMemcpyDeviceToHost));
	int count=0;
	for(int i=841659; i< 841660; i++)
	{
		printf("%d %d %d %d %d\n", squares2[i].ray_index[0], squares2[i].ray_index[1], squares2[i].ray_index[2], squares2[i].ray_index[3], squares2[i].ray_index[4]);
	}*/
	
	dim3 threadSize(512, 1, 1);
	dim3 blockSize(RayTubeNum / 512 + 1, 1, 1);
	CreateRayTubeInfo << <blockSize, threadSize >> >(d_squares, d_rays, RayTubeNum, d_rays2, d_squares_pred, d_DivRayTubeNum, lmd, fi, direction, RayTubeNumNew);
	
	/*Direction* rays2;
	rays2 = (Direction*)malloc(8000000 * sizeof(Direction));
	HANDLE_ERROR(cudaMemcpy(rays2, d_rays2, 8000000 * sizeof(Direction), cudaMemcpyDeviceToHost));
	printf("rays %d: %f %f %f\n", 841659, rays2[841659].p[0], rays2[841659].p[1], rays2[841659].p[2]);
	printf("rays %d: %f %f %f\n", 1893734, rays2[1893734].p[0], rays2[1893734].p[1], rays2[1893734].p[2]);
	printf("rays %d: %f %f %f\n", 841658, rays2[841658].p[0], rays2[841658].p[1], rays2[841658].p[2]);
	printf("rays %d: %f %f %f\n", 1683319, rays2[1683319].p[0], rays2[1683319].p[1], rays2[1683319].p[2]);
	printf("rays %d: %f %f %f\n", 2735394, rays2[2735394].p[0], rays2[2735394].p[1], rays2[2735394].p[2]);
	printf("-------------------------------------------------------\n");
	printf("rays %d: %f %f %f\n", 0, rays2[0].p[0], rays2[0].p[1], rays2[0].p[2]);
	printf("rays %d: %f %f %f\n", 1, rays2[1].p[0], rays2[1].p[1], rays2[1].p[2]);
	printf("rays %d: %f %f %f\n", 2, rays2[2].p[0], rays2[2].p[1], rays2[2].p[2]);
	printf("rays %d: %f %f %f\n", 3, rays2[3].p[0], rays2[3].p[1], rays2[3].p[2]);
	printf("rays %d: %f %f %f\n", 4, rays2[4].p[0], rays2[4].p[1], rays2[4].p[2]);
	printf("-------------------------------------------------------\n");*/
	
	DivRayTubeCenter << <SizeOfblock, SizeOfthread >> >(d_squares2, d_rays2, RayTubeNumNew, DivRayTubeNum, direction);
	
	/*//Direction* rays2;
	//rays2 = (Direction*)malloc(8000000 * sizeof(Direction));
	HANDLE_ERROR(cudaMemcpy(rays2, d_rays2, 8000000 * sizeof(Direction), cudaMemcpyDeviceToHost));
	printf("rays %d: %f %f %f\n", 841659, rays2[841659].p[0], rays2[841659].p[1], rays2[841659].p[2]);
	printf("rays %d: %f %f %f\n", 1893734, rays2[1893734].p[0], rays2[1893734].p[1], rays2[1893734].p[2]);
	printf("rays %d: %f %f %f\n", 841658, rays2[841658].p[0], rays2[841658].p[1], rays2[841658].p[2]);
	printf("rays %d: %f %f %f\n", 1683319, rays2[1683319].p[0], rays2[1683319].p[1], rays2[1683319].p[2]);
	printf("rays %d: %f %f %f\n", 2735394, rays2[2735394].p[0], rays2[2735394].p[1], rays2[2735394].p[2]);
	printf("-------------------------------------------------------\n");*/

	//CopyRayTube << <blockSize, threadSize >> >(d_squares, RayTubeNum, d_squares2, d_rays2, d_rays, d_squares_pred, d_DivRayTubeNum, RayTubeNumNew);
	
	/*//Direction* rays2;
	//rays2 = (Direction*)malloc(8000000 * sizeof(Direction));
	HANDLE_ERROR(cudaMemcpy(rays2, d_rays2, 8000000 * sizeof(Direction), cudaMemcpyDeviceToHost));
	printf("rays %d: %f %f %f\n", 841659, rays2[841659].p[0], rays2[841659].p[1], rays2[841659].p[2]);
	printf("rays %d: %f %f %f\n", 1893734, rays2[1893734].p[0], rays2[1893734].p[1], rays2[1893734].p[2]);
	printf("rays %d: %f %f %f\n", 841658, rays2[841658].p[0], rays2[841658].p[1], rays2[841658].p[2]);
	printf("rays %d: %f %f %f\n", 1683319, rays2[1683319].p[0], rays2[1683319].p[1], rays2[1683319].p[2]);
	printf("rays %d: %f %f %f\n", 2735394, rays2[2735394].p[0], rays2[2735394].p[1], rays2[2735394].p[2]);
	printf("-------------------------------------------------------\n");
	printf("rays %d: %f %f %f\n", 841659, rays2[841659].dir[0], rays2[841659].dir[1], rays2[841659].dir[2]);
	printf("rays %d: %f %f %f\n", 1893734, rays2[1893734].dir[0], rays2[1893734].dir[1], rays2[1893734].dir[2]);
	printf("rays %d: %f %f %f\n", 841658, rays2[841658].dir[0], rays2[841658].dir[1], rays2[841658].dir[2]);
	printf("rays %d: %f %f %f\n", 1683319, rays2[1683319].dir[0], rays2[1683319].dir[1], rays2[1683319].dir[2]);
	printf("rays %d: %f %f %f\n", 2735394, rays2[2735394].dir[0], rays2[2735394].dir[1], rays2[2735394].dir[2]);
	printf("-------------------------------------------------------\n");*/

}