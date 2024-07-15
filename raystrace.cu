#include "common_struct.h"
#include <math.h>
#include <stdio.h>
#include "handlerror.h"
#include "raystrace.h"

#define TPB 128
__device__ int count1[181] = { 0 };
__device__ int count2[181] = { 0 };
__device__ int count3[181] = { 0 };

__device__
bool RayIntersectBox_In(KD_Node_V *node, Direction *ray, float *intersect_point)
{
	float ptOnPlane[3];                //射线与包围盒某一个面的交点
	float *minPoint = node->box.bmin;  //包围盒的最小点
	float *maxPoint = node->box.bmax;  //包围盒的最大点

	float3 origin = ray->p;
	float3 dir = ray->dir;

	float t;
	if (dir.x != 0.f)
	{
		if (dir.x > 0)
			t = (minPoint[0] - origin.x) / dir.x;
		else
			t = (maxPoint[0] - origin.x) / dir.x;
		if (t >= 0.f)
		{
			ptOnPlane[0] = origin.x + t*dir.x;
			ptOnPlane[1] = origin.y + t*dir.y;
			ptOnPlane[2] = origin.z + t*dir.z;
			if (minPoint[1] <= ptOnPlane[1] && ptOnPlane[1] <= maxPoint[1] && minPoint[2] <= ptOnPlane[2] && ptOnPlane[2] <= maxPoint[2])
			{
				intersect_point[0] = ptOnPlane[0];
				intersect_point[1] = ptOnPlane[1];
				intersect_point[2] = ptOnPlane[2];
				return true;
			}
		}
	}
	if (dir.y != 0.f)
	{
		if (dir.y > 0)
			t = (minPoint[1] - origin.y) / dir.y;
		else
			t = (maxPoint[1] - origin.y) / dir.y;
		if (t >= 0.f)
		{
			ptOnPlane[0] = origin.x + t*dir.x;
			ptOnPlane[1] = origin.y + t*dir.y;
			ptOnPlane[2] = origin.z + t*dir.z;
			if (minPoint[2] <= ptOnPlane[2] && ptOnPlane[2] <= maxPoint[2] && minPoint[0] <= ptOnPlane[0] && ptOnPlane[0] <= maxPoint[0])
			{
				intersect_point[0] = ptOnPlane[0];
				intersect_point[1] = ptOnPlane[1];
				intersect_point[2] = ptOnPlane[2];
				return true;
			}
		}

	}
	if (dir.z != 0.f)
	{
		if (dir.z > 0)
			t = (minPoint[2] - origin.z) / dir.z;
		else
			t = (maxPoint[2] - origin.z) / dir.z;
		if (t >= 0.f)
		{
			ptOnPlane[0] = origin.x + t*dir.x;
			ptOnPlane[1] = origin.y + t*dir.y;
			ptOnPlane[2] = origin.z + t*dir.z;
			if (minPoint[1] <= ptOnPlane[1] && ptOnPlane[1] <= maxPoint[1] && minPoint[0] <= ptOnPlane[0] && ptOnPlane[0] <= maxPoint[0])
			{
				intersect_point[0] = ptOnPlane[0];
				intersect_point[1] = ptOnPlane[1];
				intersect_point[2] = ptOnPlane[2];
				return true;
			}
		}
	}
	return false;
}

__device__
bool RayIntersectBox_Out(KD_Node_V *node, Direction *ray, float *intersect_point, int *whichface)
{
	int front_or_behind = -1;
	float ptOnPlane[3];                   //射线与包围盒某一个面的交点
	float *minPoint = node->box.bmin;     //包围盒的最小点
	float *maxPoint = node->box.bmax;     //包围盒的最大点

	float3 origin = ray->p;
	float3 dir = ray->dir;

	float t;
	if (dir.x != 0.f)
	{
        //如果找出口点的话,t为远点处的t
		if (dir.x > 0)
		{
			t = (maxPoint[0] - origin.x) / dir.x;
			front_or_behind = 1;
		}
		else
		{
			t = (minPoint[0] - origin.x) / dir.x;
			front_or_behind = 0;
		}
		if (t >= 0.f)
		{
			ptOnPlane[0] = origin.x + t*dir.x;
			ptOnPlane[1] = origin.y + t*dir.y;
			ptOnPlane[2] = origin.z + t*dir.z;
			if (minPoint[1] <= ptOnPlane[1] && ptOnPlane[1] <= maxPoint[1] && minPoint[2] <= ptOnPlane[2] && ptOnPlane[2] <= maxPoint[2])
			{
				intersect_point[0] = ptOnPlane[0];
				intersect_point[1] = ptOnPlane[1];
				intersect_point[2] = ptOnPlane[2];
				*whichface = 0 + front_or_behind;
				return true;
			}
		}
	}
	if (dir.y != 0.f)
	{
		if (dir.y > 0)
		{
			t = (maxPoint[1] - origin.y) / dir.y;
			front_or_behind = 1;
		}
		else
		{
			t = (minPoint[1] - origin.y) / dir.y;
			front_or_behind = 0;
		}
		if (t >= 0.f)
		{
			ptOnPlane[0] = origin.x + t*dir.x;
			ptOnPlane[1] = origin.y + t*dir.y;
			ptOnPlane[2] = origin.z + t*dir.z;
			if (minPoint[2] <= ptOnPlane[2] && ptOnPlane[2] <= maxPoint[2] && minPoint[0] <= ptOnPlane[0] && ptOnPlane[0] <= maxPoint[0])
			{
				intersect_point[0] = ptOnPlane[0];
				intersect_point[1] = ptOnPlane[1];
				intersect_point[2] = ptOnPlane[2];
				*whichface = 2 + front_or_behind;
				return true;
			}
		}

	}
	if (dir.z != 0.f)
	{
		if (dir.z > 0)
		{
			t = (maxPoint[2] - origin.z) / dir.z;
			front_or_behind = 1;
		}
		else
		{
			t = (minPoint[2] - origin.z) / dir.z;
			front_or_behind = 0;
		}
		if (t >= 0.f)
		{
			ptOnPlane[0] = origin.x + t*dir.x;
			ptOnPlane[1] = origin.y + t*dir.y;
			ptOnPlane[2] = origin.z + t*dir.z;
			if (minPoint[1] <= ptOnPlane[1] && ptOnPlane[1] <= maxPoint[1] && minPoint[0] <= ptOnPlane[0] && ptOnPlane[0] <= maxPoint[0])
			{
				intersect_point[0] = ptOnPlane[0];
				intersect_point[1] = ptOnPlane[1];
				intersect_point[2] = ptOnPlane[2];
				*whichface = 4 + front_or_behind;
				return true;
			}
		}
	}
	return false;
}

__device__
bool RayIntersectTriangle(Direction *ray, float *P0, float *P1, float *P2, float *distance, float *intersect_point)
{
	//三角形法线  N=(P1-P0)x(P2-P0)
	float N[3];
	N[0] = (P1[1] - P0[1])*(P2[2] - P0[2]) - (P2[1] - P0[1])*(P1[2] - P0[2]);
	N[1] = (P1[2] - P0[2])*(P2[0] - P0[0]) - (P2[2] - P0[2])*(P1[0] - P0[0]);
	N[2] = (P1[0] - P0[0])*(P2[1] - P0[1]) - (P2[0] - P0[0])*(P1[1] - P0[1]);

	//三角形平面到原点距离 ，平面方程 N*P+d=0;
	float d = (-N[0] * P0[0]) + (-N[1] * P0[1]) + (-N[2] * P0[2]);
	if (N[0] * ray->dir.x + N[1] * ray->dir.y + N[2] * ray->dir.z == 0)
	{
		return false;
	}

	//射线方程 O+t*dir=P(t);
	float t = (-d - (N[0] * ray->p.x + N[1] * ray->p.y + N[2] * ray->p.z)) /
		(N[0] * ray->dir.x + N[1] * ray->dir.y + N[2] * ray->dir.z);
	if (t <= 0)
		return false;

	//射线与三角形平面交点
	float P_intersect[3];
	P_intersect[0] = ray->p.x + t*ray->dir.x;
	P_intersect[1] = ray->p.y + t*ray->dir.y;
	P_intersect[2] = ray->p.z + t*ray->dir.z;

	//计算交点是否在三角形内，参考博客https://blog.csdn.net/ZJU_fish1996/article/details/52276987
	float R[3];
	R[0] = P_intersect[0] - P0[0]; R[1] = P_intersect[1] - P0[1]; R[2] = P_intersect[2] - P0[2];
	float Q1[3];
	Q1[0] = P1[0] - P0[0]; Q1[1] = P1[1] - P0[1]; Q1[2] = P1[2] - P0[2];
	float Q2[3];
	Q2[0] = P2[0] - P0[0]; Q2[1] = P2[1] - P0[1]; Q2[2] = P2[2] - P0[2];

	float fm = (Q1[0] * Q1[0] + Q1[1] * Q1[1] + Q1[2] * Q1[2])*(Q2[0] * Q2[0] + Q2[1] * Q2[1] + Q2[2] * Q2[2]) -
		(Q1[0] * Q2[0] + Q1[1] * Q2[1] + Q1[2] * Q2[2])*(Q1[0] * Q2[0] + Q1[1] * Q2[1] + Q1[2] * Q2[2]);

	float w1, w2;
	w1 = ((Q2[0] * Q2[0] + Q2[1] * Q2[1] + Q2[2] * Q2[2])*(R[0] * Q1[0] + R[1] * Q1[1] + R[2] * Q1[2]) -
		(Q1[0] * Q2[0] + Q1[1] * Q2[1] + Q1[2] * Q2[2])*(R[0] * Q2[0] + R[1] * Q2[1] + R[2] * Q2[2])) / fm;
	w2 = ((Q1[0] * Q1[0] + Q1[1] * Q1[1] + Q1[2] * Q1[2])*(R[0] * Q2[0] + R[1] * Q2[1] + R[2] * Q2[2]) -
		(Q1[0] * Q2[0] + Q1[1] * Q2[1] + Q1[2] * Q2[2])*(R[0] * Q1[0] + R[1] * Q1[1] + R[2] * Q1[2])) / fm;
        
    //备注 由于float精度的问题，会造成部分求交失败
	if (w1 >= 0 && w2 >= 0 && w1 + w2 <= 1)
	{

		intersect_point[0] = P_intersect[0];
		intersect_point[1] = P_intersect[1];
		intersect_point[2] = P_intersect[2];
		*distance = sqrt((ray->p.x - intersect_point[0])*(ray->p.x - intersect_point[0])
			+ (ray->p.y - intersect_point[1])*(ray->p.y - intersect_point[1])
			+ (ray->p.z - intersect_point[2])*(ray->p.z - intersect_point[2]));
	}
	else
	{
		return  false;
	}
	return true;
}

__device__
void RayReflectAndUpdateRayDir(Direction *ray, float *P0, float *P1, float *P2)
{
	//求法线
	float N[3];
	N[0] = (P1[1] - P0[1])*(P2[2] - P0[2]) - (P2[1] - P0[1])*(P1[2] - P0[2]);
	N[1] = (P1[2] - P0[2])*(P2[0] - P0[0]) - (P2[2] - P0[2])*(P1[0] - P0[0]);
	N[2] = (P1[0] - P0[0])*(P2[1] - P0[1]) - (P2[0] - P0[0])*(P1[1] - P0[1]);

	//法线归一化
	float NLen = sqrt(pow(N[0], 2) + pow(N[1], 2) + pow(N[2], 2));
	N[0] = N[0] / NLen;
	N[1] = N[1] / NLen;
	N[2] = N[2] / NLen;
	if (N[0] * ray->dir.x + N[1] * ray->dir.y + N[2] * ray->dir.z <= 0)
	{
		N[0] = -N[0]; N[1] = -N[1]; N[2] = -N[2];
	}
	float LxN = N[0] * ray->dir.x + N[1] * ray->dir.y + N[2] * ray->dir.z;
	ray->dir.x = ray->dir.x - 2 * LxN * N[0];
	ray->dir.y = ray->dir.y - 2 * LxN * N[1];
	ray->dir.z = ray->dir.z - 2 * LxN * N[2];
}



__device__
void SingleRayTrace(KD_Node_V *root, Direction *ray, Prim_Box *arrays, Element *points, Triangle *triangles, float water_line)
{
	float intersect_point[3];
    //与根节点求交失败
	if (false == RayIntersectBox_In(root, ray, intersect_point))
	{
		return;
	}
	KD_Node_V leaf_node = root[0];
	int max_cycle_count=1000;          //设置最大循环次数，防止未知原因造成死循环
	while (max_cycle_count--)
	{
        //循环遍历KDTree 直到找到底层包围盒
		while ((leaf_node.IsLeaf | leaf_node.IsEmpty)==false)
		{
			if ((intersect_point[leaf_node.Split_Axis] - leaf_node.SplitPos.point[leaf_node.Split_Axis])>0.000001)
			{
				leaf_node = root[leaf_node.RightIndex];
			}
			else /*if ((intersect_point[leaf_node.Split_Axis] - leaf_node.SplitPos.point[leaf_node.Split_Axis])<0.000001)*/
			{
				leaf_node = root[leaf_node.LeftIndex];
			}
		}
		
		//到达底层包围盒
		if (leaf_node.IsEmpty)
		{
			int whichface = -1;
            
			if (false == RayIntersectBox_Out(&leaf_node, ray, intersect_point, &whichface))
			{                     
                //从包围盒出去的时候出错了，理论不该有，防止由于精度等位置原因造成这种错误。
				return;
			}
			int out_index = leaf_node.RopeIndex[whichface];
			if (out_index != -1)
			{
                //空节点，穿过
				leaf_node = root[leaf_node.RopeIndex[whichface]];
			}
			else
			{
				//到达出口，退出
				return;
			}
		}
		else
		{
			//进入非空叶节点
			int calced_triangle_id = -1;
			while (ray->times <= 2)
			{
				float near_triangle_distance = 999999999;
				int intersect_triangle_id = -1;
				float current_intersect_point[3];
				for (int triIndex = leaf_node.begin; triIndex <= leaf_node.end; triIndex++)
				{
					if (triIndex == calced_triangle_id)
						continue;
					float P0[3], P1[3], P2[3];
					float intersect_point_on_triangle[3];
					float distance;
					P0[0] = points[triangles[arrays[triIndex].Box_Index].Points[0]].point[0];
					P1[0] = points[triangles[arrays[triIndex].Box_Index].Points[1]].point[0];
					P2[0] = points[triangles[arrays[triIndex].Box_Index].Points[2]].point[0];
					P0[1] = points[triangles[arrays[triIndex].Box_Index].Points[0]].point[1];
					P1[1] = points[triangles[arrays[triIndex].Box_Index].Points[1]].point[1];
					P2[1] = points[triangles[arrays[triIndex].Box_Index].Points[2]].point[1];
					P0[2] = points[triangles[arrays[triIndex].Box_Index].Points[0]].point[2];
					P1[2] = points[triangles[arrays[triIndex].Box_Index].Points[1]].point[2];
					P2[2] = points[triangles[arrays[triIndex].Box_Index].Points[2]].point[2];

					if (true == RayIntersectTriangle(ray, P0, P1, P2, &distance, intersect_point_on_triangle))
					{
                        //判断交点是否在包围盒内，如果不在，说明它在隔壁的包围盒，跳过
						if (intersect_point_on_triangle[0] >= leaf_node.box.bmin[0]
							&& intersect_point_on_triangle[0] <= leaf_node.box.bmax[0]   
							&& intersect_point_on_triangle[1] >= leaf_node.box.bmin[1]
							&& intersect_point_on_triangle[1] <= leaf_node.box.bmax[1]
							&& intersect_point_on_triangle[2] >= leaf_node.box.bmin[2]
							&& intersect_point_on_triangle[2] <= leaf_node.box.bmax[2])
						{
							if (0.01<distance && distance < near_triangle_distance)
							{
								near_triangle_distance = distance;
								intersect_triangle_id = triIndex;
								current_intersect_point[0] = intersect_point_on_triangle[0];
								current_intersect_point[1] = intersect_point_on_triangle[1];
								current_intersect_point[2] = intersect_point_on_triangle[2];
							}
						}
					}
				}
                //不等于-1的时候就说明射线与包围盒中的三角形有交点，然后就需要更新ray的原点，方向和据远点距离
				if (intersect_triangle_id != -1)
				{
					float P0[3], P1[3], P2[3];
					P0[0] = points[triangles[arrays[intersect_triangle_id].Box_Index].Points[0]].point[0];
					P1[0] = points[triangles[arrays[intersect_triangle_id].Box_Index].Points[1]].point[0];
					P2[0] = points[triangles[arrays[intersect_triangle_id].Box_Index].Points[2]].point[0];
					P0[1] = points[triangles[arrays[intersect_triangle_id].Box_Index].Points[0]].point[1];
					P1[1] = points[triangles[arrays[intersect_triangle_id].Box_Index].Points[1]].point[1];
					P2[1] = points[triangles[arrays[intersect_triangle_id].Box_Index].Points[2]].point[1];
					P0[2] = points[triangles[arrays[intersect_triangle_id].Box_Index].Points[0]].point[2];
					P1[2] = points[triangles[arrays[intersect_triangle_id].Box_Index].Points[1]].point[2];
					P2[2] = points[triangles[arrays[intersect_triangle_id].Box_Index].Points[2]].point[2];
					
					//若相交点在水线上，则抛弃之 22.11.7 jzy
					if(current_intersect_point[2]>water_line)
						break;
					
					RayReflectAndUpdateRayDir(ray, P0, P1, P2);
					ray->p.x = current_intersect_point[0];
					ray->p.y = current_intersect_point[1];
					ray->p.z = current_intersect_point[2];
					ray->distance += near_triangle_distance;
					ray->times++;
					ray->triangle_index = arrays[intersect_triangle_id].Box_Index;
					//ray->triangle_index = intersect_triangle_id;

					calced_triangle_id = intersect_triangle_id;
				}
				else
				{
                    //搜索完毕，与非空包围盒中的三角形都没有交点，退出非空包围盒
					break;
				}
			}
			if (ray->times == 3)
			{
			    //射线反射了三次，退出
				return;
			}
			int whichface;
			if (false == RayIntersectBox_Out(&leaf_node, ray, intersect_point, &whichface))
			{
                //求出口包围盒失败，如从这里退出，则是异常状况
				return;
			}
			int next_index = leaf_node.RopeIndex[whichface];
			if (next_index < 0)
			{
			     //从包围盒出来之后就出去了，正常返回
				return;
			}
				
			leaf_node = root[next_index];
		}
	}
	return;

}


//__device__ 
//bool IfValid(int a, int b, int c, int d)
//{
//	/*if ((a == b && b == c &&a != -1) || (a == b&&b == d&&a != -1) || (b == c&&c == d&&c != -1) || (a == c&&c == d&&a != -1))*/
//	if (a == b && b == c &&a != -1 && c == d)
//		return true;
//	else
//		return false;
//}
//
//
//__global__ 
//void sideraytracekernel_1d_gpu(KD_Node_V *d_root, Direction *d_rays, Square *d_squares, Prim_Box *d_array, Element *d_points, Triangle *d_triangles, int width, int height)
//{
//	int i = blockIdx.x*blockDim.x + threadIdx.x;
//	if (i >= (width + 1)*(height + 1))
//		return;
//	Direction ray = d_rays[i];
//	SingleRayTrace(d_root, &ray, d_array, d_points, d_triangles);
//	d_rays[i] = ray;
//}
//
//__global__ 
//void centraytracekernel_1d_gpu(KD_Node_V *d_root, Direction *d_rays, Square *d_squares, Prim_Box *d_array, Element *d_points, Triangle *d_triangles, int width, int height)
//{
//	int i = blockIdx.x*blockDim.x + threadIdx.x;
//
//	if (i >= width*height)
//		return;
//
//	int a = d_rays[d_squares[i].ray_index[0]].triangle_index;
//	int b = d_rays[d_squares[i].ray_index[1]].triangle_index;
//	int c = d_rays[d_squares[i].ray_index[2]].triangle_index;
//	int d = d_rays[d_squares[i].ray_index[3]].triangle_index;
//
//	bool value = IfValid(a, b, c, d);
//	if (value)
//	{
//		Direction ray = d_rays[d_squares[i].ray_index[4]];
//		SingleRayTrace(d_root, &ray, d_array, d_points, d_triangles);
//		d_rays[d_squares[i].ray_index[4]] = ray;
//	}
//	d_squares[i].right = value;
//	
//}
//
//
//void allraystrace_v2(Direction *d_rays, Square *d_squares, int width, int height,
//	KD_Node_V *d_root, Prim_Box *d_array, int prim_boxnum,
//	Element *d_points, int pointsnum, Triangle *d_triangles, int trianglesnum)
//{
//	const int gridSize1 = ((width + 1)*(height + 1) + TPB - 1) / TPB;
//	const int blockSize1 = TPB;
//
//	sideraytracekernel_1d_gpu << <gridSize1, blockSize1 >> >(d_root, d_rays, d_squares, d_array, d_points, d_triangles, width, height);
//
//	HANDLE_ERROR(cudaGetLastError());
//	//HANDLE_ERROR(cudaDeviceSynchronize());
//
//	const int gridSize2 = (width * height + TPB - 1) / TPB;
//	const int blockSize2 = TPB;
//
//	centraytracekernel_1d_gpu << <gridSize2, blockSize2>> >(d_root, d_rays, d_squares, d_array, d_points, d_triangles, width, height);
//
//	HANDLE_ERROR(cudaGetLastError());
//	//HANDLE_ERROR(cudaDeviceSynchronize());
//}

//声线管束有效性判断,四个顶点在同一个三角面元内即为有效
__device__ bool ValidRayTube(int a, int b, int c, int d)
{
	if (a == b && b == c && a != -1 && c == d)
		return true;
	else
		return false;
}

__device__ bool InvalidRayTube(int a, int b, int c, int d)
{
	return (a == b && b == c && a == -1 && c == d);
}



__global__ void centraytracekernel_1d_gpu(KD_Node_V *d_root, Direction *d_rays, Square *d_squares, Prim_Box *d_array, Element *d_points, Triangle *d_triangles,
	int DivRayTubeNum, int RayTubeNum, Axis direction, float water_line)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;

	if (i >= DivRayTubeNum)
		return;

	Direction tmp = d_rays[i + RayTubeNum];
	//tmp.caled = false;
	//tmp.flag = false;
	tmp.distance = 0;
	tmp.times = 0;
	tmp.triangle_index = -1;
	tmp.dir.x = -direction.x;
	tmp.dir.y = -direction.y;
	tmp.dir.z = -direction.z;
	SingleRayTrace(d_root, &tmp, d_array, d_points, d_triangles, water_line);
	d_rays[i + RayTubeNum] = tmp;

}

__global__ void CreateCenterRay_Div(Direction *d_rays, Square *d_squares, int RayTubeNum, int DivRayTubeNum, int* d_DivRayTubeNum, int* d_squares_pred)
{
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	int tid = threadIdx.x;

	if (idx >= RayTubeNum)
		return;

	if (d_squares[idx].IsDivRayTube == true)
	{
		int index = d_DivRayTubeNum[blockIdx.x] + d_squares_pred[idx] + RayTubeNum * 2 + DivRayTubeNum;
		d_squares[idx].CenterRayIndex = index;
		d_rays[index].p = d_squares[idx].CenterRay;
		//d_rays[index].p.x = d_squares[idx].CenterRay.x;
		//d_rays[index].p.y = d_squares[idx].CenterRay.y;
		//d_rays[index].p.z = d_squares[idx].CenterRay.z;
	}
}

__global__ void CreateCenterRay(Direction *d_rays, Square *d_squares, int RayTubeNum, int* d_DivRayTubeNum, int* d_squares_pred)
{
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	int tid = threadIdx.x;

	if (idx >= RayTubeNum)
		return;

	if (d_squares[idx].IsDivRayTube == true)
	{
		int index = d_DivRayTubeNum[blockIdx.x] + d_squares_pred[idx] + RayTubeNum;
		d_squares[idx].CenterRayIndex = index;
		d_rays[index].p = d_squares[idx].CenterRay;
		//d_rays[index].p[0] = d_squares[idx].CenterRay.x;
		//d_rays[index].p[1] = d_squares[idx].CenterRay.y;
		//d_rays[index].p[2] = d_squares[idx].CenterRay.z;
	}
}

__global__ void ExclusiveSumScan(int* d_squares_pred, int RayTubeNum)
{
	int tid = threadIdx.x;
	int idx = 2 * tid + 2 * blockDim.x * blockIdx.x;

	int offset = 1;
	int num = (RayTubeNum + 2 * blockDim.x - 1) / (2 * blockDim.x); //长度为array_length的数组被blockDim.x*2分割成了num个segment
	int len = num * 2 * blockDim.x;//补齐后数组长度为len

	__shared__ int temp[512];

	//边界填充
	for (int i = tid; i < blockDim.x; i += blockDim.x)
	{
		temp[2 * tid] = 0;
		temp[2 * tid + 1] = 0;
	}
	__syncthreads();

	if (idx < RayTubeNum)
	{
		temp[2 * tid] = d_squares_pred[idx];
	}

	if (idx + 1 < RayTubeNum)
	{
		temp[2 * tid + 1] = d_squares_pred[idx + 1];
	}
	__syncthreads();

	//up-sweep phase上行阶段
	for (int j = blockDim.x; j > 0; j >>= 1)
	{
		if (tid < j)
		{
			int ai = offset * (2 * tid + 1) - 1;
			int bi = offset * (2 * tid + 2) - 1;
			temp[bi] += temp[ai];
		}
		offset *= 2;
		__syncthreads();
	}

	//down-sweep phase下行阶段
	if (tid == 0)
		temp[511] = 0;

	for (int j = 1; j < (blockDim.x * 2); j *= 2)
	{
		offset >>= 1;
		__syncthreads();
		if (tid < j)
		{
			int ai = offset * (2 * tid + 1) - 1;
			int bi = offset * (2 * tid + 2) - 1;

			float t = temp[ai];
			temp[ai] = temp[bi];
			temp[bi] += t;
		}
		__syncthreads();
	}

	if (idx < RayTubeNum)
	{
		d_squares_pred[idx] = temp[2 * tid];
	}
	if (idx + 1 < RayTubeNum)
	{
		d_squares_pred[idx + 1] = temp[2 * tid + 1];
	}
}

__global__ void ScanInArray(Square *d_squares, int* d_squares_pred, int RayTubeNum)
{
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	int tid = threadIdx.x;

	__shared__ bool sDivFlag[512];
	sDivFlag[tid] = (idx < RayTubeNum) ? d_squares[idx].IsDivRayTube : false;
	__syncthreads();

	if (idx < RayTubeNum)
	{
		d_squares_pred[idx] = (sDivFlag[tid] == true) ? 1 : 0;
	}
}

__global__ void BlellochScan1(int* d_DivRayTubeNum, int array_length, int* d_sum_gmem, int* d_sum_Gmem, int* count1_s)
{
	int tid = threadIdx.x;
	int idx = 2 * tid + 2 * blockDim.x * blockIdx.x;
	int id = tid + blockDim.x * blockIdx.x;
	
	int offset = 1;
	int num = (array_length + 2 * blockDim.x - 1) / (2 * blockDim.x); //长度为array_length的数组被blockDim.x*2分割成了num个segment
	int len = num * 2 * blockDim.x;//补齐后数组长度为len

	__shared__ int temp[512];

	//边界填充
	for (int i = tid; i < blockDim.x; i += blockDim.x)
	{
		temp[2 * tid] = 0;
		temp[2 * tid + 1] = 0;
	}
	__syncthreads();

	if (idx < array_length)
	{
		temp[2 * tid] = d_DivRayTubeNum[idx];
	}

	if (idx + 1 < array_length)
	{
		temp[2 * tid + 1] = d_DivRayTubeNum[idx + 1];
	}
	__syncthreads();

	//up-sweep phase上行阶段
	for (int j = blockDim.x; j > 0; j >>= 1)
	{
		if (tid < j)
		{
			int ai = offset * (2 * tid + 1) - 1;
			int bi = offset * (2 * tid + 2) - 1;
			temp[bi] += temp[ai];
		}
		offset *= 2;
		__syncthreads();
	}

	//down-sweep phase下行阶段
	__shared__ int islast;
	if (tid == 0)
	{
		d_sum_gmem[blockIdx.x + 1] = temp[2 * blockDim.x - 1];
		temp[2 * blockDim.x - 1] = 0;
		__threadfence();

		int value = atomicAdd(count1_s, 1);
		islast = (value == gridDim.x - 1);
	}
	__syncthreads();

	if (islast)
	{
		int i = id - blockIdx.x * blockDim.x;

		if (i < gridDim.x)
		{
			int tmp = 0;
			for (int n = 0; n < i + 2; n++)
			{
				tmp += d_sum_gmem[n];
			}
			//__syncthreads();
			d_sum_Gmem[i + 1] = tmp;
		}
		__syncthreads();
	}

	for (int j = 1; j < (blockDim.x * 2); j *= 2)
	{
		offset >>= 1;
		__syncthreads();
		if (tid < j)
		{
			int ai = offset * (2 * tid + 1) - 1;
			int bi = offset * (2 * tid + 2) - 1;

			float t = temp[ai];
			temp[ai] = temp[bi];
			temp[bi] += t;
		}
		__syncthreads();
	}

	if (idx < array_length)
	{
		//d_DivRayTubeNum[idx] = temp[2 * tid] + d_sum_gmem[blockIdx.x];
		d_DivRayTubeNum[idx] = temp[2 * tid];
	}
	if (idx + 1 < array_length)
	{
		//d_DivRayTubeNum[idx + 1] = temp[2 * tid + 1] + d_sum_gmem[blockIdx.x];
		d_DivRayTubeNum[idx + 1] = temp[2 * tid + 1];
	}
}

__global__ void DivRayTubeNumAdd(int* d_DivRayTubeNum, int* d_sum_gmem, int array_length)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < array_length)
	{
		d_DivRayTubeNum[idx] += d_sum_gmem[blockIdx.x];
	}
}

__global__ void IsValidRayTube_gpu(Direction *d_rays, Square *d_squares, int RayTubeNum, int* d_DivRayTubeNum)
{

	__shared__ bool sDivFlag[512];

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int tid = threadIdx.x;

	if (idx >= RayTubeNum)
		return;

	int a = d_rays[d_squares[idx].CornerRayIndex.x].triangle_index;
	int b = d_rays[d_squares[idx].CornerRayIndex.y].triangle_index;
	int c = d_rays[d_squares[idx].CornerRayIndex.z].triangle_index;
	int d = d_rays[d_squares[idx].CornerRayIndex.w].triangle_index;
	d_squares[idx].right = ValidRayTube(a, b, c, d);//参与积分的声线管束
	d_squares[idx].IsDivRayTube = !(InvalidRayTube(a, b, c, d) ^ d_squares[idx].right);//需要分裂的声线管束

	sDivFlag[tid] = (idx < RayTubeNum) ? d_squares[idx].IsDivRayTube : false;
	__syncthreads();

	if (sDivFlag[tid] == true)
	{
		atomicAdd(&d_DivRayTubeNum[blockIdx.x], 1);
	}
}

__global__ void sideraytracekernel_1d_gpu(KD_Node_V *d_root, Direction *d_rays, Square *d_squares, Prim_Box *d_array, Element *d_points, Triangle *d_triangles,
	int totalraysnum, Axis direction, float water_line)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= totalraysnum)
		return;

	Direction tmp = d_rays[i];//2019-04-08
	//tmp.caled = false;
	//tmp.flag = false;
	tmp.distance = 0;
	tmp.times = 0;
	tmp.triangle_index = -1;
	tmp.dir.x = -direction.x;
	tmp.dir.y = -direction.y;
	tmp.dir.z = -direction.z;
	SingleRayTrace(d_root, &tmp, d_array, d_points, d_triangles, water_line);
	d_rays[i] = tmp;
}

void allraystrace_v2(Direction *d_rays, Square *d_squares, int width, int height, KD_Node_V *d_root, Prim_Box *d_array, Element *d_points, Triangle *d_triangles,
	int* d_DivRayTubeNum, int* DivRayTubeNum, int* d_sum_gmem, int* d_sum_Gmem, int* d_squares_pred, Axis direction, float angle, float water_line)
{
	const int gridSize1 = ((width + 1) * (height + 1) + TPB - 1) / TPB;
	const int blockSize1 = TPB;

	float temp = sqrt((direction.x * direction.x) + (direction.y * direction.y) + (direction.z * direction.z));
	direction.x /= temp;
	direction.y /= temp;
	direction.z /= temp;

	sideraytracekernel_1d_gpu << <gridSize1, blockSize1 >> >(d_root, d_rays, d_squares, d_array, d_points, d_triangles, (width + 1) * (height + 1), direction, water_line);
	HANDLE_ERROR(cudaGetLastError());
	HANDLE_ERROR(cudaDeviceSynchronize());

	dim3 threadSize(512, 1, 1);
	dim3 blockSize(width * height / 512 + 1, 1, 1);
	int array_length = blockSize.x;

	IsValidRayTube_gpu << <blockSize, threadSize >> >(d_rays, d_squares, width * height, d_DivRayTubeNum);
	HANDLE_ERROR(cudaGetLastError());
	//HANDLE_ERROR(cudaDeviceSynchronize());

	dim3 ThreadSize(256, 1, 1);
	dim3 BlockSize(array_length / 512 + 1, 1, 1);
	int ang = (int)angle;
	int x[38575];
	cudaMemcpy(x,d_DivRayTubeNum,sizeof(int)*38575,cudaMemcpyDeviceToHost);
	
	int *count1_s;
	cudaMalloc((void**)&count1_s,sizeof(int));
	int value = 0;
	cudaMemcpy(count1_s,&value,sizeof(int),cudaMemcpyHostToDevice);
	BlellochScan1 << <BlockSize, ThreadSize >> >(d_DivRayTubeNum, array_length, d_sum_gmem, d_sum_Gmem, count1_s);
	HANDLE_ERROR(cudaGetLastError());
	cudaFree(count1_s);
	
	//HANDLE_ERROR(cudaDeviceSynchronize());

/*
	int* sum_gmem;
	sum_gmem = (int*)malloc(35 * sizeof(int));
	printf("sum_gmem!\n");
	HANDLE_ERROR(cudaMemcpy(sum_gmem, d_sum_gmem, 35 * sizeof(int), cudaMemcpyDeviceToHost));
	for(int i=0; i<35; i++)
	{
	printf("%d ",sum_gmem[i]);
	}
	printf("\n");
	free(sum_gmem);
*/

	int DivRayTubeNumIdx = BlockSize.x;
	HANDLE_ERROR(cudaMemcpy(DivRayTubeNum, &d_sum_Gmem[DivRayTubeNumIdx], sizeof(int), cudaMemcpyDeviceToHost));
	//debug*********
	//printf("DivRayTubeNumIdx=%d\n",DivRayTubeNumIdx);
	//***********
	DivRayTubeNumAdd << <BlockSize, threadSize >> >(d_DivRayTubeNum, d_sum_Gmem, array_length);
	HANDLE_ERROR(cudaGetLastError());
	HANDLE_ERROR(cudaDeviceSynchronize());

	ScanInArray << <blockSize, threadSize >> >(d_squares, d_squares_pred, width * height);
	HANDLE_ERROR(cudaGetLastError());
	
	dim3 Sizethread(256, 1, 1);
	dim3 Sizeblock(width * height / 512 + 1, 1, 1);
	ExclusiveSumScan << <Sizeblock, Sizethread >> >(d_squares_pred, width * height);
	HANDLE_ERROR(cudaGetLastError());

	CreateCenterRay << <blockSize, threadSize >> >(d_rays, d_squares, (width + 1) * (height + 1), d_DivRayTubeNum, d_squares_pred);
	HANDLE_ERROR(cudaGetLastError());
	const int gridSize2 = (*DivRayTubeNum + TPB - 1) / TPB;
	const int blockSize2 = TPB;
	centraytracekernel_1d_gpu << <gridSize2, blockSize2 >> >(d_root, d_rays, d_squares, d_array, d_points, d_triangles, *DivRayTubeNum, (width + 1) * (height + 1), direction, water_line);
	//HANDLE_ERROR(cudaGetLastError());	
	HANDLE_ERROR(cudaDeviceSynchronize());

}


__global__ void BlellochScan2(int* d_DivRayTubeNum, int array_length, int* d_sum_gmem, int* d_sum_Gmem, int* count2_s)
{
	int tid = threadIdx.x;
	int idx = 2 * tid + 2 * blockDim.x * blockIdx.x;
	int id = tid + blockDim.x * blockIdx.x;

	int offset = 1;
	int num = (array_length + 2 * blockDim.x - 1) / (2 * blockDim.x); //长度为array_length的数组被blockDim.x*2分割成了num个segment
	int len = num * 2 * blockDim.x;//补齐后数组长度为len

	__shared__ int temp[512];

	//边界填充
	for (int i = tid; i < blockDim.x; i += blockDim.x)
	{
		temp[2 * tid] = 0;
		temp[2 * tid + 1] = 0;
	}
	__syncthreads();

	if (idx < array_length)
	{
		temp[2 * tid] = d_DivRayTubeNum[idx];
	}

	if (idx + 1 < array_length)
	{
		temp[2 * tid + 1] = d_DivRayTubeNum[idx + 1];
	}
	__syncthreads();

	//up-sweep phase上行阶段
	for (int j = blockDim.x; j > 0; j >>= 1)
	{
		if (tid < j)
		{
			int ai = offset * (2 * tid + 1) - 1;
			int bi = offset * (2 * tid + 2) - 1;
			temp[bi] += temp[ai];
		}
		offset *= 2;
		__syncthreads();
	}

	//down-sweep phase下行阶段
	__shared__ int islast;
	if (tid == 0)
	{
		d_sum_gmem[blockIdx.x + 1] = temp[2 * blockDim.x - 1];
		temp[2 * blockDim.x - 1] = 0;
		__threadfence();

		int value = atomicAdd(count2_s, 1);
		islast = (value == gridDim.x - 1);
	}
	__syncthreads();

	if (islast)
	{
		int i = id - blockIdx.x * blockDim.x;
		if (i < gridDim.x)
		{
			int tmp = 0;
			for (int n = 0; n < i + 2; n++)
			{
				tmp += d_sum_gmem[n];
			}
			//__syncthreads();
			d_sum_Gmem[i + 1] = tmp;
		}
		__syncthreads();
	}

	for (int j = 1; j < (blockDim.x * 2); j *= 2)
	{
		offset >>= 1;
		__syncthreads();
		if (tid < j)
		{
			int ai = offset * (2 * tid + 1) - 1;
			int bi = offset * (2 * tid + 2) - 1;

			float t = temp[ai];
			temp[ai] = temp[bi];
			temp[bi] += t;
		}
		__syncthreads();
	}

	if (idx < array_length)
	{
		d_DivRayTubeNum[idx] = temp[2 * tid];
	}
	if (idx + 1 < array_length)
	{
		d_DivRayTubeNum[idx + 1] = temp[2 * tid + 1];
	}
}

__global__ void sideraytracekernel_DivRay_gpu(KD_Node_V *d_root, Direction *d_rays, Prim_Box *d_array, Element *d_points, Triangle *d_triangles,
	int d_totalDivRayNum, Axis direction, float water_line)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i >= d_totalDivRayNum)
		return;

	Direction tmp = d_rays[i];//2019-04-08
	//tmp.caled = false;
	//tmp.flag = false;
	tmp.distance = 0;
	tmp.times = 0;
	tmp.triangle_index = -1;
	tmp.dir.x = -direction.x;
	tmp.dir.y = -direction.y;
	tmp.dir.z = -direction.z;
	SingleRayTrace(d_root, &tmp, d_array, d_points, d_triangles, water_line);
	d_rays[i] = tmp;
}

__global__ void centraytracekernel_DivRay_gpu(KD_Node_V *d_root, Direction *d_rays, Square *d_squares, Prim_Box *d_array, Element *d_points, Triangle *d_triangles,
	int DivRayTubeNum, int d_totalDivRayNum, Axis direction, float water_line)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i >= DivRayTubeNum)
		return;

	Direction temp = d_rays[i + d_totalDivRayNum * 2 + DivRayTubeNum];
	//temp.caled = false;
	//temp.flag = false;
	temp.distance = 0;
	temp.times = 0;
	temp.triangle_index = -1;
	temp.dir.x = -direction.x;
	temp.dir.y = -direction.y;
	temp.dir.z = -direction.z;
	SingleRayTrace(d_root, &temp, d_array, d_points, d_triangles, water_line);
	d_rays[i + d_totalDivRayNum * 2 + DivRayTubeNum] = temp;
}

void allraystrace_DivRayFirst(Direction *d_rays2, Square *d_squares2, int d_totalDivRayNum, KD_Node_V *d_root, Prim_Box *d_array, Element *d_points, Triangle *d_triangles,
	int* d_DivRayTubeNum, int* DivRayTubeNum, int* d_sum_gmem, int* d_sum_Gmem, int* d_squares_pred, Axis direction, float angle, float water_line)
{
	d_totalDivRayNum *= 4;
	const int gridSize1 = (d_totalDivRayNum + TPB - 1) / TPB;
	const int blockSize1 = TPB;

	float temp = sqrt((direction.x * direction.x) + (direction.y * direction.y) + (direction.z * direction.z));
	direction.x /= temp;
	direction.y /= temp;
	direction.z /= temp;

	sideraytracekernel_DivRay_gpu << <gridSize1, blockSize1 >> >(d_root, d_rays2, d_array, d_points, d_triangles, d_totalDivRayNum, direction, water_line);
	HANDLE_ERROR(cudaGetLastError());
	//HANDLE_ERROR(cudaDeviceSynchronize());


	dim3 threadSize(512, 1, 1);
	dim3 blockSize(d_totalDivRayNum / 512 + 1, 1, 1);
	int array_length = blockSize.x;

	IsValidRayTube_gpu << <blockSize, threadSize >> >(d_rays2, d_squares2, d_totalDivRayNum, d_DivRayTubeNum);

	dim3 ThreadSize(256, 1, 1);
	dim3 BlockSize(array_length / 512 + 1, 1, 1);
	int ang = (int)angle;
	
	int *count2_s;
	cudaMalloc((void**)&count2_s,sizeof(int));
	int value = 0;
	cudaMemcpy(count2_s,&value,sizeof(int),cudaMemcpyHostToDevice);
	BlellochScan2 << <BlockSize, ThreadSize >> >(d_DivRayTubeNum, array_length, d_sum_gmem, d_sum_Gmem, count2_s);
	cudaFree(count2_s);

	int DivRayTubeNumIdx = BlockSize.x;
	cudaMemcpy(DivRayTubeNum, &d_sum_Gmem[DivRayTubeNumIdx], sizeof(int), cudaMemcpyDeviceToHost);

	DivRayTubeNumAdd << <BlockSize, threadSize >> >(d_DivRayTubeNum, d_sum_Gmem, array_length);

	//int DivRayTubeNumIdx = BlockSize.x;

	ScanInArray << <blockSize, threadSize >> >(d_squares2, d_squares_pred, d_totalDivRayNum);

	dim3 Sizethread(256, 1, 1);
	dim3 Sizeblock(d_totalDivRayNum / 512 + 1, 1, 1);
	ExclusiveSumScan << <Sizeblock, Sizethread >> >(d_squares_pred, d_totalDivRayNum);

	//cudaMemcpy(DivRayTubeNum, &d_sum_gmem[DivRayTubeNumIdx], sizeof(int), cudaMemcpyDeviceToHost);

	CreateCenterRay_Div << <blockSize, threadSize >> >(d_rays2, d_squares2, d_totalDivRayNum, *DivRayTubeNum, d_DivRayTubeNum, d_squares_pred);

	const int gridSize2 = (*DivRayTubeNum + TPB - 1) / TPB;
	const int blockSize2 = TPB;

	centraytracekernel_DivRay_gpu << <gridSize2, blockSize2 >> >(d_root, d_rays2, d_squares2, d_array, d_points, d_triangles, *DivRayTubeNum, d_totalDivRayNum, direction, water_line);

	HANDLE_ERROR(cudaGetLastError());
	//HANDLE_ERROR(cudaDeviceSynchronize());
}



__global__ void BlellochScan3(int* d_DivRayTubeNum, int array_length, int* d_sum_gmem, int* d_sum_Gmem, int* count3_s)
{
	int tid = threadIdx.x;
	int idx = 2 * tid + 2 * blockDim.x * blockIdx.x;
	int id = tid + blockDim.x * blockIdx.x;

	int offset = 1;
	int num = (array_length + 2 * blockDim.x - 1) / (2 * blockDim.x); //长度为array_length的数组被blockDim.x*2分割成了num个segment
	int len = num * 2 * blockDim.x;//补齐后数组长度为len

	__shared__ int temp[512];

	//边界填充
	for (int i = tid; i < blockDim.x; i += blockDim.x)
	{
		temp[2 * tid] = 0;
		temp[2 * tid + 1] = 0;
	}
	__syncthreads();

	if (idx < array_length)
	{
		temp[2 * tid] = d_DivRayTubeNum[idx];
	}

	if (idx + 1 < array_length)
	{
		temp[2 * tid + 1] = d_DivRayTubeNum[idx + 1];
	}
	__syncthreads();

	//up-sweep phase上行阶段
	for (int j = blockDim.x; j > 0; j >>= 1)
	{
		if (tid < j)
		{
			int ai = offset * (2 * tid + 1) - 1;
			int bi = offset * (2 * tid + 2) - 1;
			temp[bi] += temp[ai];
		}
		offset *= 2;
		__syncthreads();
	}

	//down-sweep phase下行阶段
	__shared__ int islast;
	if (tid == 0)
	{
		d_sum_gmem[blockIdx.x + 1] = temp[2 * blockDim.x - 1];
		temp[2 * blockDim.x - 1] = 0;
		__threadfence();

		int value = atomicAdd(count3_s, 1);
		islast = (value == gridDim.x - 1);
	}
	__syncthreads();

	if (islast)
	{
		int i = id - blockIdx.x * blockDim.x;
		if (i < gridDim.x)
		{
			int tmp = 0;
			for (int n = 0; n < i + 2; n++)
			{
				tmp += d_sum_gmem[n];
			}
			//__syncthreads();
			d_sum_Gmem[i + 1] = tmp;
		}
		__syncthreads();
	}

	for (int j = 1; j < (blockDim.x * 2); j *= 2)
	{
		offset >>= 1;
		__syncthreads();
		if (tid < j)
		{
			int ai = offset * (2 * tid + 1) - 1;
			int bi = offset * (2 * tid + 2) - 1;

			float t = temp[ai];
			temp[ai] = temp[bi];
			temp[bi] += t;
		}
		__syncthreads();
	}

	if (idx < array_length)
	{
		d_DivRayTubeNum[idx] = temp[2 * tid];
	}
	if (idx + 1 < array_length)
	{
		d_DivRayTubeNum[idx + 1] = temp[2 * tid + 1];
	}
}

void allraystrace_DivRaySecond(Direction *d_rays2, Square *d_squares2, int d_totalDivRayNum, KD_Node_V *d_root, Prim_Box *d_array, Element *d_points, Triangle *d_triangles,
	int* d_DivRayTubeNum, int* DivRayTubeNum, int* d_sum_gmem, int* d_sum_Gmem, int* d_squares_pred, Axis direction, float angle, float water_line)
{
	d_totalDivRayNum *= 4;
	const int gridSize1 = (d_totalDivRayNum + TPB - 1) / TPB;
	const int blockSize1 = TPB;

	float temp = sqrt((direction.x * direction.x) + (direction.y * direction.y) + (direction.z * direction.z));
	direction.x /= temp;
	direction.y /= temp;
	direction.z /= temp;

	sideraytracekernel_DivRay_gpu << <gridSize1, blockSize1 >> >(d_root, d_rays2, d_array, d_points, d_triangles, d_totalDivRayNum, direction, water_line);
	HANDLE_ERROR(cudaGetLastError());
	//HANDLE_ERROR(cudaDeviceSynchronize());


	dim3 threadSize(512, 1, 1);
	dim3 blockSize(d_totalDivRayNum / 512 + 1, 1, 1);
	int array_length = blockSize.x;

	IsValidRayTube_gpu << <blockSize, threadSize >> >(d_rays2, d_squares2, d_totalDivRayNum, d_DivRayTubeNum);

	dim3 ThreadSize(256, 1, 1);
	dim3 BlockSize(array_length / 512 + 1, 1, 1);
	int ang = (int)angle;
	int *count3_s;
	cudaMalloc((void**)&count3_s,sizeof(int));
	int value = 0;
	cudaMemcpy(count3_s,&value,sizeof(int),cudaMemcpyHostToDevice);
	BlellochScan3 << <BlockSize, ThreadSize >> >(d_DivRayTubeNum, array_length, d_sum_gmem, d_sum_Gmem, count3_s);
	cudaFree(count3_s);

	int DivRayTubeNumIdx = BlockSize.x;
	cudaMemcpy(DivRayTubeNum, &d_sum_Gmem[DivRayTubeNumIdx], sizeof(int), cudaMemcpyDeviceToHost);

	DivRayTubeNumAdd << <BlockSize, threadSize >> >(d_DivRayTubeNum, d_sum_Gmem, array_length);

	//int DivRayTubeNumIdx = BlockSize.x;

	ScanInArray << <blockSize, threadSize >> >(d_squares2, d_squares_pred, d_totalDivRayNum);

	dim3 Sizethread(256, 1, 1);
	dim3 Sizeblock(d_totalDivRayNum / 512 + 1, 1, 1);
	ExclusiveSumScan << <Sizeblock, Sizethread >> >(d_squares_pred, d_totalDivRayNum);

	//cudaMemcpy(DivRayTubeNum, &d_sum_gmem[DivRayTubeNumIdx], sizeof(int), cudaMemcpyDeviceToHost);

	CreateCenterRay_Div << <blockSize, threadSize >> >(d_rays2, d_squares2, d_totalDivRayNum, *DivRayTubeNum, d_DivRayTubeNum, d_squares_pred);

	const int gridSize2 = (*DivRayTubeNum + TPB - 1) / TPB;
	const int blockSize2 = TPB;

	centraytracekernel_DivRay_gpu << <gridSize2, blockSize2 >> >(d_root, d_rays2, d_squares2, d_array, d_points, d_triangles, *DivRayTubeNum, d_totalDivRayNum, direction, water_line);

	HANDLE_ERROR(cudaGetLastError());
	//HANDLE_ERROR(cudaDeviceSynchronize());
}

__global__ void centraytracekernel_DivRaySecond_gpu(Direction *d_rays, Square *d_squares, int d_totalDivRayNum)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;

	if (i >= d_totalDivRayNum)
		return;

	int a = d_rays[d_squares[i].CornerRayIndex.x].triangle_index;
	int b = d_rays[d_squares[i].CornerRayIndex.y].triangle_index;
	int c = d_rays[d_squares[i].CornerRayIndex.z].triangle_index;
	int d = d_rays[d_squares[i].CornerRayIndex.w].triangle_index;
	d_squares[i].right = ValidRayTube(a, b, c, d);

}

void allraystrace_DivRayThird(Direction *d_rays2, Square *d_squares2, int d_totalDivRayNum, KD_Node_V *d_root, Prim_Box *d_array,
	Element *d_points, Triangle *d_triangles, Axis direction, float water_line)
{
	d_totalDivRayNum *= 4;

	/*Square* squares2;
	squares2 = (Square*)malloc(15976009 * sizeof(Square));
	HANDLE_ERROR(cudaMemcpy(squares2, d_squares2, 15976009 * sizeof(Square), cudaMemcpyDeviceToHost));
	int count=0;
	for(int i=9716; i< 9717; i++)
	{
	printf("%d %d %d %d %d\n", squares2[i].ray_index[0], squares2[i].ray_index[1], squares2[i].ray_index[2], squares2[i].ray_index[3], squares2[i].ray_index[4]);
	}
	Direction* rays2;
	rays2 = (Direction*)malloc(31960013 * sizeof(Direction));
	HANDLE_ERROR(cudaMemcpy(rays2, d_rays2, 31960013 * sizeof(Direction), cudaMemcpyDeviceToHost));
	printf("rays %d: %.7f %.7f %.7f\n", 9574360, rays2[9574360].p[0], rays2[9574360].p[1], rays2[9574360].p[2]);
	printf("rays %d: %.7f %.7f %.7f\n", 9716, rays2[9716].p[0], rays2[9716].p[1], rays2[9716].p[2]);
	printf("rays %d: %.7f %.7f %.7f\n", 19131717, rays2[19131717].p[0], rays2[19131717].p[1], rays2[19131717].p[2]);
	printf("rays %d: %.7f %.7f %.7f\n", 9719, rays2[9719].p[0], rays2[9719].p[1], rays2[9719].p[2]);
	//printf("rays %d: %.7f %.7f %.7f\n", 2735394, rays2[2735394].p[0], rays2[2735394].p[1], rays2[2735394].p[2]);
	printf("-------------------------------------------------------\n");*/

	const int gridSize = (d_totalDivRayNum + TPB - 1) / TPB;
	const int blockSize = TPB;

	float temp = sqrt((direction.x * direction.x) + (direction.y * direction.y) + (direction.z * direction.z));
	direction.x /= temp;
	direction.y /= temp;
	direction.z /= temp;

	sideraytracekernel_DivRay_gpu << <gridSize, blockSize >> >(d_root, d_rays2, d_array, d_points, d_triangles, d_totalDivRayNum, direction, water_line);

	/*HANDLE_ERROR(cudaMemcpy(rays2, d_rays2, 31960013 * sizeof(Direction), cudaMemcpyDeviceToHost));
	printf("rays %d: %.7f %.7f %.7f\n", 9574360, rays2[9574360].p[0], rays2[9574360].p[1], rays2[9574360].p[2]);
	printf("rays %d: %.7f %.7f %.7f\n", 9716, rays2[9716].p[0], rays2[9716].p[1], rays2[9716].p[2]);
	printf("rays %d: %.7f %.7f %.7f\n", 19131717, rays2[19131717].p[0], rays2[19131717].p[1], rays2[19131717].p[2]);
	printf("rays %d: %.7f %.7f %.7f\n", 9719, rays2[9719].p[0], rays2[9719].p[1], rays2[9719].p[2]);
	//printf("rays %d: %.7f %.7f %.7f\n", 2735394, rays2[2735394].p[0], rays2[2735394].p[1], rays2[2735394].p[2]);
	printf("-------------------------------------------------------\n");*/

	HANDLE_ERROR(cudaGetLastError());
	//HANDLE_ERROR(cudaDeviceSynchronize());

	centraytracekernel_DivRaySecond_gpu << <gridSize, blockSize >> >(d_rays2, d_squares2, d_totalDivRayNum);

	HANDLE_ERROR(cudaGetLastError());
	//HANDLE_ERROR(cudaDeviceSynchronize());
}
