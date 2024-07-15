#include "simple_time.h"
#include "raystrace.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "device_functions.h"
#include "tree2vector.h"
#include <math.h>
#include <stdio.h>
#include "handlerror.h"
#include "virtualface_gpu.h"
#include"integral_gpu.h"

void MallocOnGPU(int dwidth, int dheight, Direction** d_rays1, Square** d_squares1, Direction** d_rays2, Square** d_squares2,
	RayBeamInfo** d_effrays, Vector** d_center, Vector** d_axis, MatStruct** d_transMat, ReimOutput** d_reim,
	float** d_sum_re, float** d_sum_im, float** d_sum_sre,float** d_sum_sim,int** d_DivRayTubeNum, int** d_sum_gmem, int** d_sum_Gmem, int** d_squares_pred)
{
	//printf("dwidth_malloc=%d, dheight_malloc=%d\n",dwidth,dheight);

	HANDLE_ERROR(cudaMalloc((void**)d_rays1, sizeof(Direction)*((dwidth + 1)*(dheight + 1) + dwidth*dheight)));//0.82GB//4850000
	HANDLE_ERROR(cudaMalloc((void**)d_rays2, sizeof(Direction)*((dwidth + 1)*(dheight + 1) + dwidth*dheight+4850000)));//0.82GB
	HANDLE_ERROR(cudaMalloc((void**)d_squares1, sizeof(Square)*dwidth*dheight));//0.46GB
	HANDLE_ERROR(cudaMalloc((void**)d_squares2, sizeof(Square)*dwidth*dheight));//0.46GB

	//printf("sizeof d_rays1:%d\n",sizeof(Direction)*((dwidth + 1)*(dheight + 1) + dwidth*dheight));\
	//printf("sizeof d_squares1:%d\n",sizeof(Square)*dwidth*dheight);

	int raysBeamNum = dwidth * dheight;
	//printf("raysBeamNum malloc = %d\n",raysBeamNum);
	//printf("raysNum malloc = %d\n",(dwidth + 1)*(dheight + 1) + dwidth*dheight);
	dim3 blockSize((raysBeamNum + 511)/ 512, 1, 1);
	
	//3.12GB
	HANDLE_ERROR(cudaMalloc((void**)d_effrays, raysBeamNum * sizeof(RayBeamInfo)));//1.64GB
	HANDLE_ERROR(cudaMalloc((void**)d_center, raysBeamNum * sizeof(Vector)));//0.11GB
	HANDLE_ERROR(cudaMalloc((void**)d_axis, 3 * raysBeamNum * sizeof(Vector)));//0.33GB
	HANDLE_ERROR(cudaMalloc((void**)d_transMat, raysBeamNum * sizeof(MatStruct)));//0.45GB
	HANDLE_ERROR(cudaMalloc((void**)d_reim, raysBeamNum * sizeof(ReimOutput)));//0.075GB
	HANDLE_ERROR(cudaMalloc((void**)d_sum_re, blockSize.x * sizeof(float)));//0.076MB
	HANDLE_ERROR(cudaMalloc((void**)d_sum_im, blockSize.x * sizeof(float)));//0.076MB
	HANDLE_ERROR(cudaMalloc((void**)d_sum_sre, blockSize.x * sizeof(float)));//0.076MB
	HANDLE_ERROR(cudaMalloc((void**)d_sum_sim, blockSize.x * sizeof(float)));//0.076MB

	HANDLE_ERROR(cudaMalloc((void**)d_DivRayTubeNum, blockSize.x * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)d_sum_gmem, 35 * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)d_sum_Gmem, 35 * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)d_squares_pred, raysBeamNum * sizeof(int)));

	/*HANDLE_ERROR(cudaMemset(*d_effrays, 0, raysBeamNum * sizeof(RayBeamInfo)));
	HANDLE_ERROR(cudaMemset(*d_center, 0, raysBeamNum * sizeof(Vector)));
	HANDLE_ERROR(cudaMemset(*d_axis, 0, 3 * raysBeamNum * sizeof(Vector)));
	HANDLE_ERROR(cudaMemset(*d_transMat, 0, raysBeamNum * sizeof(MatStruct)));
	HANDLE_ERROR(cudaMemset(*d_reim, 0, raysBeamNum * sizeof(comp)));
	HANDLE_ERROR(cudaMemset(*d_sum_re, 0, blockSize.x * sizeof(float)));
	HANDLE_ERROR(cudaMemset(*d_sum_im, 0, blockSize.x * sizeof(float)));
	HANDLE_ERROR(cudaMemset(*d_DivRayTubeNum, 0, blockSize.x * sizeof(int)));
	HANDLE_ERROR(cudaMemset(*d_sum_gmem, 0, 35 * sizeof(int)));
	HANDLE_ERROR(cudaMemset(*d_sum_Gmem, 0, 35 * sizeof(int)));
	HANDLE_ERROR(cudaMemset(*d_squares_pred, 0, raysBeamNum * sizeof(int)));*/
}

void MemsetOnGPU1(int dwidth, int dheight, Direction** d_rays1, Square** d_squares1, RayBeamInfo** d_effrays, Vector** d_center, Vector** d_axis, 
	MatStruct** d_transMat, ReimOutput** d_reim, float** d_sum_re, float** d_sum_im,
	int** d_DivRayTubeNum, int** d_sum_gmem, int** d_sum_Gmem, int** d_squares_pred)
{
	int raysBeamNum = dwidth * dheight;
	dim3 threadSize(512, 1, 1);
	dim3 blockSize(raysBeamNum / 512 + 1, 1, 1);

	HANDLE_ERROR(cudaMemset(*d_squares1, 0, dwidth*dheight * sizeof(Square)));
	HANDLE_ERROR(cudaMemset(*d_rays1, 0, ((dwidth + 1)*(dheight + 1) + dwidth*dheight) * sizeof(Direction)));
	HANDLE_ERROR(cudaMemset(*d_effrays, 0, raysBeamNum * sizeof(RayBeamInfo)));
	HANDLE_ERROR(cudaMemset(*d_center, 0, raysBeamNum * sizeof(Vector)));
	HANDLE_ERROR(cudaMemset(*d_axis, 0, 3 * raysBeamNum * sizeof(Vector)));
	HANDLE_ERROR(cudaMemset(*d_transMat, 0, raysBeamNum * sizeof(MatStruct)));
	HANDLE_ERROR(cudaMemset(*d_reim, 0, raysBeamNum * sizeof(ReimOutput)));
	HANDLE_ERROR(cudaMemset(*d_sum_re, 0, blockSize.x * sizeof(float)));
	HANDLE_ERROR(cudaMemset(*d_sum_im, 0, blockSize.x * sizeof(float)));

	HANDLE_ERROR(cudaMemset(*d_DivRayTubeNum, 0, blockSize.x * sizeof(int)));
	HANDLE_ERROR(cudaMemset(*d_sum_gmem, 0, 35 * sizeof(int)));
	HANDLE_ERROR(cudaMemset(*d_sum_Gmem, 0, 35 * sizeof(int)));
	HANDLE_ERROR(cudaMemset(*d_squares_pred, 0, raysBeamNum * sizeof(int)));
}

void MemsetOnGPU2(int dwidth, int dheight, Direction** d_rays2, Square** d_squares2, RayBeamInfo** d_effrays, Vector** d_center, Vector** d_axis, 
	MatStruct** d_transMat, ReimOutput** d_reim, float** d_sum_re, float** d_sum_im,
	int** d_DivRayTubeNum, int** d_sum_gmem, int** d_sum_Gmem, int** d_squares_pred)
{
	int raysBeamNum = dwidth * dheight;
	dim3 threadSize(512, 1, 1);
	dim3 blockSize(raysBeamNum / 512 + 1, 1, 1);

	HANDLE_ERROR(cudaMemset(*d_squares2, 0, dwidth*dheight * sizeof(Square)));
	HANDLE_ERROR(cudaMemset(*d_rays2, 0, ((dwidth + 1)*(dheight + 1) + dwidth*dheight+4850000) * sizeof(Direction)));
	HANDLE_ERROR(cudaMemset(*d_effrays, 0, raysBeamNum * sizeof(RayBeamInfo)));
	HANDLE_ERROR(cudaMemset(*d_center, 0, raysBeamNum * sizeof(Vector)));
	HANDLE_ERROR(cudaMemset(*d_axis, 0, 3 * raysBeamNum * sizeof(Vector)));
	HANDLE_ERROR(cudaMemset(*d_transMat, 0, raysBeamNum * sizeof(MatStruct)));
	HANDLE_ERROR(cudaMemset(*d_reim, 0, raysBeamNum * sizeof(ReimOutput)));
	HANDLE_ERROR(cudaMemset(*d_sum_re, 0, blockSize.x * sizeof(float)));
	HANDLE_ERROR(cudaMemset(*d_sum_im, 0, blockSize.x * sizeof(float)));
	HANDLE_ERROR(cudaMemset(*d_DivRayTubeNum, 0, blockSize.x * sizeof(int)));
	HANDLE_ERROR(cudaMemset(*d_sum_gmem, 0, 35 * sizeof(int)));
	HANDLE_ERROR(cudaMemset(*d_sum_Gmem, 0, 35 * sizeof(int)));
	HANDLE_ERROR(cudaMemset(*d_squares_pred, 0, raysBeamNum * sizeof(int)));
}

void MemsetOnGPU3(int dwidth, int dheight, RayBeamInfo** d_effrays, Vector** d_center, Vector** d_axis, MatStruct** d_transMat, ReimOutput** d_reim, float** d_sum_re, float** d_sum_im)
{
	int raysBeamNum = dwidth * dheight;
	dim3 threadSize(512, 1, 1);
	dim3 blockSize(raysBeamNum / 512 + 1, 1, 1);

	HANDLE_ERROR(cudaMemset(*d_effrays, 0, raysBeamNum * sizeof(RayBeamInfo)));
	HANDLE_ERROR(cudaMemset(*d_center, 0, raysBeamNum * sizeof(Vector)));
	HANDLE_ERROR(cudaMemset(*d_axis, 0, 3 * raysBeamNum * sizeof(Vector)));
	HANDLE_ERROR(cudaMemset(*d_transMat, 0, raysBeamNum * sizeof(MatStruct)));
	HANDLE_ERROR(cudaMemset(*d_reim, 0, raysBeamNum * sizeof(ReimOutput)));
	HANDLE_ERROR(cudaMemset(*d_sum_re, 0, blockSize.x * sizeof(float)));
	HANDLE_ERROR(cudaMemset(*d_sum_im, 0, blockSize.x * sizeof(float)));
}

void MemsetOnGPU(int dwidth, int dheight, Direction** d_rays1, Square** d_squares1, Direction** d_rays2, Square** d_squares2,
	RayBeamInfo** d_effrays, Vector** d_center, Vector** d_axis, MatStruct** d_transMat, ReimOutput** d_reim,
	float** d_sum_re, float** d_sum_im, float** d_sum_sre,float** d_sum_sim, int** d_DivRayTubeNum, int** d_sum_gmem, int** d_sum_Gmem, int** d_squares_pred)
{
	int raysBeamNum = dwidth * dheight;
	dim3 threadSize(512, 1, 1);
	dim3 blockSize(raysBeamNum / 512 + 1, 1, 1);

	HANDLE_ERROR(cudaMemset(*d_squares2, 0, dwidth*dheight * sizeof(Square)));
	HANDLE_ERROR(cudaMemset(*d_rays2, 0, ((dwidth + 1)*(dheight + 1) + dwidth*dheight+4850000) * sizeof(Direction)));
	HANDLE_ERROR(cudaMemset(*d_squares1, 0, dwidth*dheight * sizeof(Square)));
	HANDLE_ERROR(cudaMemset(*d_rays1, 0, ((dwidth + 1)*(dheight + 1) + dwidth*dheight) * sizeof(Direction)));
	HANDLE_ERROR(cudaMemset(*d_effrays, 0, raysBeamNum * sizeof(RayBeamInfo)));
	HANDLE_ERROR(cudaMemset(*d_center, 0, raysBeamNum * sizeof(Vector)));
	HANDLE_ERROR(cudaMemset(*d_axis, 0, 3 * raysBeamNum * sizeof(Vector)));
	HANDLE_ERROR(cudaMemset(*d_transMat, 0, raysBeamNum * sizeof(MatStruct)));
	HANDLE_ERROR(cudaMemset(*d_reim, 0, raysBeamNum * sizeof(ReimOutput)));
	HANDLE_ERROR(cudaMemset(*d_sum_re, 0, blockSize.x * sizeof(float)));
	HANDLE_ERROR(cudaMemset(*d_sum_im, 0, blockSize.x * sizeof(float)));
	HANDLE_ERROR(cudaMemset(*d_sum_sre, 0, blockSize.x * sizeof(float)));//20200919
HANDLE_ERROR(cudaMemset(*d_sum_sim, 0, blockSize.x * sizeof(float)));//20200919
	HANDLE_ERROR(cudaMemset(*d_DivRayTubeNum, 0, blockSize.x * sizeof(int)));
	HANDLE_ERROR(cudaMemset(*d_sum_gmem, 0, 35 * sizeof(int)));
	HANDLE_ERROR(cudaMemset(*d_sum_Gmem, 0, 35 * sizeof(int)));
	HANDLE_ERROR(cudaMemset(*d_squares_pred, 0, raysBeamNum * sizeof(int)));
}

void free_virtualface(Direction *rays, Square *squares)
{
	cudaFree(rays);
	cudaFree(squares);
}

void FreeOnGPU(Direction* d_rays1, Square* d_squares1, Direction* d_rays2, Square* d_squares2, RayBeamInfo* d_effrays, Vector* d_center,
	Vector* d_axis, MatStruct* d_transMat, ReimOutput* d_reim, float* d_sum_re, float* d_sum_im, float* d_sum_sre,float* d_sum_sim, int* d_DivRayTubeNum, int* d_sum_gmem, int* d_sum_Gmem, int* d_squares_pred)
{
	cudaFree(d_rays1);
	cudaFree(d_rays2);
	cudaFree(d_squares1);
	cudaFree(d_squares2);
	
	cudaFree(d_center);
	cudaFree(d_axis);
	cudaFree(d_transMat);
	cudaFree(d_effrays);
	cudaFree(d_reim);
	cudaFree(d_sum_re);
	cudaFree(d_sum_im);
	cudaFree(d_sum_sre);
cudaFree(d_sum_sim);
	cudaFree(d_DivRayTubeNum);
	cudaFree(d_sum_gmem);
	cudaFree(d_sum_Gmem);
	cudaFree(d_squares_pred);
}

void free_data(Prim_Box *d_array, KD_Node_V *d_root, Element *d_points, Triangle *d_triangles)
{
	cudaFree(d_array);
	cudaFree(d_root);
	cudaFree(d_points);
	cudaFree(d_triangles);
}

void Destroy_Tree(KD_Node *root)
{
	if (root)
	{
		Destroy_Tree(root->LeftChild);
		Destroy_Tree(root->RightChild);
	}
	free(root);
	root = NULL;
}
