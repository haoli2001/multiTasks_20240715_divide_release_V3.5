#include<cuda_runtime.h>
#include<device_launch_parameters.h>
#include<device_functions.h>
#include<string.h>
#include "scalfuc.h"
#include "handlerror.h"
#include "sm_20_atomic_functions.h"
#define  D2R 3.14159265358979/180.0
__global__ void reduce_add_sre(float* d_sum_re, RayBeamInfo* rays, int raysBeamNum)
{
	__shared__ float sdata[512];

	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	int tid = threadIdx.x;

	sdata[tid] = (idx < raysBeamNum) ? rays[idx].parameter.p[1] : 0;
	__syncthreads();

	for (int s = blockDim.x / 2; s > 0; s >>= 1)
	{
		if (tid < s) {
			sdata[tid] += sdata[tid + s];
		}
		__syncthreads();
	}

	if (tid == 0)
	{
		d_sum_re[blockIdx.x] = sdata[0];
	}
}
__global__ void reduce_sre(float* d_sum_re, float* d_in, int Num)
{
	__shared__ float sdata[512];

	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	int tid = threadIdx.x;

	sdata[tid] = (idx < Num) ? d_in[idx] : 0;
	__syncthreads();

	for (int s = blockDim.x / 2; s > 0; s >>= 1)
	{
		if (tid < s) {
			sdata[tid] += sdata[tid + s];
		}
		__syncthreads();
	}

	if (tid == 0)
	{
		d_sum_re[blockIdx.x] = sdata[0];
	}
}
__global__ void reduce_add_sim(float* d_sum_re, RayBeamInfo* rays, int raysBeamNum)
{
	__shared__ float sdata[512];

	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	int tid = threadIdx.x;

	sdata[tid] = (idx < raysBeamNum) ? rays[idx].parameter.p[2] : 0;
	__syncthreads();

	for (int s = blockDim.x / 2; s > 0; s >>= 1)
	{
		if (tid < s) {
			sdata[tid] += sdata[tid + s];
		}
		__syncthreads();
	}

	if (tid == 0)
	{
		d_sum_re[blockIdx.x] = sdata[0];
	}
}
__global__ void reduce_sim(float* d_sum_re, float* d_in, int Num)
{
	__shared__ float sdata[512];

	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	int tid = threadIdx.x;

	sdata[tid] = (idx < Num) ? d_in[idx] : 0;
	__syncthreads();

	for (int s = blockDim.x / 2; s > 0; s >>= 1)
	{
		if (tid < s) {
			sdata[tid] += sdata[tid + s];
		}
		__syncthreads();
	}

	if (tid == 0)
	{
		d_sum_re[blockIdx.x] = sdata[0];
	}
}


//以下四个函数cdiv cmul cabs d_ReflectCoeff_2为反射系数计算需要  姬梓遇20210913
#define A2R 3.1415926535897/180

__device__ __host__ comp cdiv(comp z1,comp z2)       
{ double x1,x2,y1,y2;
  comp z;
  x1=z1.re;
  x2=z2.re;
  y1=z1.im;
  y2=z2.im;
  z.re=(x1*x2+y1*y2)/(x2*x2+y2*y2);
  z.im=(x2*y1-y2*x1)/(x2*x2+y2*y2);
  return z;
}

__device__ __host__ comp cmul(comp z1,comp z2)
{ double x1,x2,y1,y2;
  comp z;
  x1=z1.re;
  x2=z2.re;
  y1=z1.im;
  y2=z2.im;
  z.re=x1*x2-y1*y2;
  z.im=x1*y2+y1*x2;
  return z;
}

__device__ __host__ double cabs(comp z)                      
{
	double x,y;
	x=z.re;
	y=z.im;
	return sqrt(x*x+y*y);
}

//计算双层层敷瓦T体的反射系数
__device__ float d_ReflectCoeff_2(float f, float theta)
{
	double theta0, theta1, theta2;
	double c0, c1, c2, c3;
	double rou0, rou1, rou2, rou3;
	double Ee2;
	double Z0, Z3, Z4;
	double eta1, eta2;
	double k1, k2, k3;
	double d1, d2, d3;
	double h1, h2, h3;
	double R, AR;
	double phi;
	double AR0, AR1;
	comp temp1, temp2, Z23, Z12, Z01, Zin, cw1, cw2, Z1, Z2;

	theta0 = theta;  //入射角度
	rou0 = 1000;     //水密度（单位kg/m^3）
	c0 = 1500;       //水中声速（单位m/s）
	Z0 = rou0 * c0;  //水的阻抗

	rou1 = 1039;     //橡胶（介质1）密度
	c1 = 1470;     //介质1中的等效声速
	eta1 = 0.4;      //损耗因子
	phi = atan(eta1);
	cw1.re = c1 * pow((1 * 1 + eta1*eta1), 0.25)*cos(phi / 2);
	cw1.im = c1 * pow((1 * 1 + eta1*eta1), 0.25)*sin(phi / 2);//粘弹材料波速
	Z1.re = rou1 * cw1.re;    //材料的纵波波阻抗
	Z1.im = rou1 * cw1.im;    //材料的纵波波阻抗
	k1 = 2 * 180 * A2R * f / c1;//介质1波数

	rou2 = 1090;     //橡胶（介质1）密度
	eta2 = 0.5;      //损耗因子
	phi = atan(eta2);
	Ee2 = 1e9;      //介质2的杨氏模量
	c2 = sqrt(Ee2 / rou2);//介质2中的等效声速
	cw2.re = c2 * pow((1 * 1 + eta2*eta2), 0.25)*cos(phi / 2);
	cw2.im = c2 * pow((1 * 1 + eta2*eta2), 0.25)*sin(phi / 2);//粘弹材料波速
	Z2.re = rou2 * cw2.re;    //材料的纵波波阻抗
	Z2.im = rou2 * cw2.im;    //材料的纵波波阻抗
	k2 = 2 * 180 * A2R * f / c2;//介质波数

	rou3 = 7850;     //钢（衬底）密度
	c3 = 5200;       //介质2中的等效声速
	Z3 = rou3 * c3;  //衬底(钢的阻抗大于水的20倍)
	k3 = 2 * 180 * A2R * f / c3;//介质波数,衬底为钢板时的传播速度

	d1 = 2e-3;       //介质1厚度
	d2 = 2e-3;       //介质2厚度
	d3 = 3e-3;       //衬底厚度

					 //////根据Snell折射定理计算各个介质层中的入射角度
	AR0 = sin(theta0) * c0 / c1;  //定义全反射系数
	if (fabs(AR0) > 1)
	{
		theta1 = 90 * A2R;
	}
	else
	{
		theta1 = asin(sin(theta0) * c0 / c1);
	}

	AR1 = sin(theta1) * c1 / c2;  //定义全反射系数
	if (fabs(AR1) > 1)
	{
		theta2 = 90 * A2R;
	}
	else
	{
		theta2 = asin(sin(theta1) * c1 / c2);
	}
	h1 = tan(k1 * cos(theta0) * d1);
	h2 = tan(k2 * cos(theta1) * d2);
	h3 = tan(k3 * cos(theta2) * d3);
	Z4 = Z0;

	temp1.im = Z3 * h3;
	temp1.re = Z4;
	temp2.im = Z4 * h3;
	temp2.re = Z3;
	Z23 = cdiv(temp1, temp2);
	Z23.re = Z3 * Z23.re;
	Z23.im = Z3 * Z23.im;


	temp1.im = Z23.im + Z2.re * h2;
	temp1.re = Z23.re - h2*Z2.im;
	temp2.im = Z23.re * h2 + Z2.im;
	temp2.re = Z2.re - Z23.im * h2;
	Z12 = cdiv(temp1, temp2);
	Z12 = cmul(Z12, Z2);

	temp1.im = Z12.im + Z1.re * h1;
	temp1.re = Z12.re - h1*Z1.im;
	temp2.im = Z1.im + Z12.re * h1;
	temp2.re = Z1.re - Z12.im * h1;
	Z01 = cdiv(temp1, temp2);
	Z01 = cmul(Z01, Z1);
	Zin = Z01;
	//Zin = Z1 * (Z3 * (Z2 - h1 * h2 * Z1) + jay * Z2 * (Z2 * h2 + h1 * Z1))/(Z2 * (Z1 - h1 * h2 * Z2) + jay * Z3 * (Z1 * h2 + h1 * Z2));  //输入阻抗
	R = fabs((cabs(Zin) - Z0) / (cabs(Zin) + Z0)); //反射系数
	return R;
}

__global__ void transMat(RayBeamInfo* rays, int raysBeamNum,  int ig, ConfigStruct config, Vector* d_center, Axis_slx New_receive_points) {
	float CSpeed = 1500.0;//shengsu
	float fend = config.time_end_frequency;//姬梓遇20210831
	float fbeg = config.time_start_frequency;//姬梓遇20210831
	float Tao = config.tao/1000.0; //单位:s 姬梓遇20210831
	float fs = config.sampling_frequency;//caiyanglv 姬梓遇20210831

	float velocity1 = config.velocity1;//目标速度1 姬梓遇
	float velocity2 = config.velocity2;//目标速度1 姬梓遇
    float velocity12 = velocity1 - velocity2;//相对投影速度 LV 20220720
	float dopGene = (CSpeed + velocity12) / (CSpeed - velocity12);//多普勒压缩因子
    float band = fend - fbeg;
	if(fabs(band) > 0.001)  //宽带信号则考虑脉冲宽度压缩
		Tao = Tao / dopGene;
	fbeg = config.time_start_frequency * dopGene;//调制后起始频率
	fend = fbeg + band * dopGene;
	band = fend - fbeg;
	float K0 = band / Tao;//wangying	
	//float Z0 ;//询问意义ZO->intel.cu
	int taosize = int(fs * Tao);
	//float wavenum=4.188;
	float start_alpha;

	if(config.continue_alpha == -1)
		start_alpha = config.start_alpha;
	else
		start_alpha = config.continue_alpha;
	
	
	//float tao0= ig  * (1 / fs);
	float wavenum= 2.0 * PI * fbeg / CSpeed;
    float f_ig= 0;
	float reflect_coeff= 1;
	//PointCoor point0, point1, point2;
	//Vector normal_vector;
	int pendZeroNum = 0;//wangying
	comp integralConst1 = { 0,0 };
	comp integralConst2 = { 0,0 };
	comp integralConst = { 0,0 };

	//float preTime = (ig - 1) * (1 / fs);
	//float furTime= (ig + 1) * (1 / fs);
    
    float tao0 = 0;  //snw
	float preTime =0 ;
	float furTime= 0;


	float Vector_proj_mod;
    float alpha00 = 0.0;
	for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < raysBeamNum; idx += gridDim.x * blockDim.x)//2018-12-11
	{

		if (rays[idx].p_cent_distance == 0)//20190116
			return;
        
		//*********************对每一块面片投影多普勒进行计算***********************20221229 snw//
		/*
        Axis_slx Vector_vobj;
		Vector_vobj.p[0] = 1.0 * config.velocity2 / cos(start_alpha * D2R - PI);
		Vector_vobj.p[1] = 0.0 * config.velocity2 / cos(start_alpha * D2R - PI);
		Vector_vobj.p[2] = 0.0 * config.velocity2 / cos(start_alpha * D2R - PI);

		Axis_slx Vector_proj;
	
		Vector_proj.p[0] = d_center[idx].x - New_receive_points.p[0];
		Vector_proj.p[1] = d_center[idx].y - New_receive_points.p[1];
		Vector_proj.p[2] = d_center[idx].z - New_receive_points.p[2];
		Vector_proj_mod = sqrt(Vector_proj.p[0]*Vector_proj.p[0] + Vector_proj.p[1]*Vector_proj.p[1] + Vector_proj.p[2]*Vector_proj.p[2]);
		Vector_proj.p[0] = Vector_proj.p[0] / Vector_proj_mod;
		Vector_proj.p[1] = Vector_proj.p[1] / Vector_proj_mod;
		Vector_proj.p[2] = Vector_proj.p[2] / Vector_proj_mod;
		// 修改前 LV 20230107 
		Axis_slx Vector_vtor;
		alpha00 = PI - start_alpha * D2R - config.H_axis_angle * D2R;     // 修改后 LV 20230707 
		//alpha00 = - (acos(Vector_proj.p[0]) - config.H_axis_angle * D2R); 
		Vector_vtor.p[0] = cos(alpha00) * config.velocity1 / cos(config.H_axis_angle * D2R);
		Vector_vtor.p[1] = 0.0;
		Vector_vtor.p[2] = sin(alpha00) * config.velocity1 / cos(config.H_axis_angle * D2R);

		// ****************
		velocity1 = Vector_vtor.p[0]*Vector_proj.p[0] + Vector_vtor.p[1]*Vector_proj.p[1] + Vector_vtor.p[2]*Vector_proj.p[2];
		velocity2 = Vector_vobj.p[0]*Vector_proj.p[0] + Vector_vobj.p[1]*Vector_proj.p[1] + Vector_vobj.p[2]*Vector_proj.p[2];

		velocity12 = velocity1 - velocity2;//相对投影速度 LV 20220720
		dopGene = (CSpeed + velocity12) / (CSpeed - velocity12);//多普勒压缩因子
		band = config.time_end_frequency - config.time_start_frequency;
		if(fabs(band) > 0.001)              //宽带信号则考虑脉冲宽度压缩
		{
			Tao = config.tao / 1000.0 / dopGene;    // 20221230 错误Tao 替换config.tao/1000
		}

		fbeg = config.time_start_frequency * dopGene;//调制后起始频率
		fend = fbeg + (config.time_end_frequency - config.time_start_frequency) * dopGene; // LV20230106 错误 * dopGene * dopGene ， 应该只有一个 * dopGene
		band = fend - fbeg;
		K0 = band / Tao;

		taosize = int(fs * Tao);*/

		//*********************对每一块面片投影多普勒进行计算***********************// 


		float recv_p_cent_distance = sqrt(pow(d_center[idx].x - New_receive_points.p[0], 2) + pow(d_center[idx].y - New_receive_points.p[1], 2) + pow(d_center[idx].z - New_receive_points.p[2], 2));
		pendZeroNum = int( (recv_p_cent_distance + rays[idx].p_cent_distance) / CSpeed * fs) ;//wangying  snw:pendZeroNum从2倍距离算起 jzy:发射点距离+接收点距离	
		//pendZeroNum = int(2 * rays[idx].p_cent_distance / CSpeed * fs) ;//wangying  snw:pendZeroNum从2倍距离算起
        
       f_ig = fbeg + (ig - pendZeroNum)  * K0 * (1.0 / fs);
       tao0= (ig - pendZeroNum) * (1.0 / fs);//snw
	   preTime = (ig - pendZeroNum - 1) * (1.0 / fs);
	   furTime= (ig - pendZeroNum + 1) * (1.0 / fs);


		if ((ig<pendZeroNum) || (ig>(pendZeroNum + taosize))) //wangying
		{
			rays[idx].parameter.p[0] = 0;
			rays[idx].parameter.p[1] = 0.0;
			rays[idx].parameter.p[2] = 0.0;
		} 
		else {
		/*	Z0 = 1.0f / 4.0f * (rays[idx].ray_index[0].p[0] + rays[idx].ray_index[1].p[0] + rays[idx].ray_index[2].p[0] + rays[idx].ray_index[3].p[0]);
			point0 = rays[idx].ray_index[0];
			point1 = rays[idx].ray_index[1];
			point2 = rays[idx].ray_index[2];
			normal_vector.x = (point0.p[1] - point2.p[1]) * (point1.p[2] - point0.p[2]) - (point1.p[1] - point0.p[1]) * (point0.p[2] - point2.p[2]);
			normal_vector.y = (point0.p[2] - point2.p[2]) * (point1.p[0] - point0.p[0]) - (point1.p[2] - point0.p[2]) * (point0.p[0] - point2.p[0]);
			normal_vector.z = (point0.p[0] - point2.p[0]) * (point1.p[1] - point0.p[1]) - (point1.p[0] - point0.p[0]) * (point0.p[1] - point2.p[1]);
			float gorden = sqrt(normal_vector.x * normal_vector.x + normal_vector.x * normal_vector.x + normal_vector.x * normal_vector.x);*///  gorden->intel.cu

			//rays[idx].parameter.p[0] = gorden;*/
            // LV 20221230 发现 wavenum 计算问题 需修改（原有代码）	
			wavenum = 2 * PI * fbeg / CSpeed;
			//wavenum= 2 * PI * f_ig / CSpeed;// LV 20221230 发现 wavenum 计算问题 修改（开始）
			integralConst1.re = cos(-2.0 * PI * (fbeg * preTime + 0.5 * K0 * tao0 * preTime) + 2.0 * wavenum * rays[idx].Z0); // 20221230 引入固定直流分量？？？  2 * 2 * pi * f_ig / CSpeed * rays[idx].Z0 ??
			integralConst1.im = sin(-2.0 * PI * (fbeg * preTime + 0.5 * K0 * tao0 * preTime) + 2.0 * wavenum * rays[idx].Z0); 

			integralConst2.re = cos(-2.0 * PI * (fbeg * furTime + 0.5 * K0 * tao0 * furTime) + 2.0 * wavenum * rays[idx].Z0);
			integralConst2.im = sin(-2.0 * PI * (fbeg * furTime + 0.5 * K0 * tao0 * furTime) + 2.0 * wavenum * rays[idx].Z0);

			/*wavenum = 2 * PI * f_ig / CSpeed;
			//wavenum= 2 * PI * f_ig / CSpeed;// LV 20221230 发现 wavenum 计算问题 修改（开始）
			integralConst1.re = cos(-2 * PI * (fbeg * preTime + 0.5 * K0 * ig / fs * preTime) + 2 * wavenum * rays[idx].Z0); // 20221230 引入固定直流分量？？？  2 * 2 * pi * f_ig / CSpeed * rays[idx].Z0 ??
			integralConst1.im = sin(-2 * PI * (fbeg * preTime + 0.5 * K0 * ig / fs * preTime) + 2 * wavenum * rays[idx].Z0); 

			integralConst2.re = cos(-2 * PI * (fbeg * furTime + 0.5 * K0 * ig / fs * furTime) + 2 * wavenum * rays[idx].Z0);
			integralConst2.im = sin(-2 * PI * (fbeg * furTime + 0.5 * K0 * ig / fs * furTime) + 2 * wavenum * rays[idx].Z0);*/
            
			//integralConst.re = (integralConst1.re - integralConst2.re) * rays[idx].gorden / rays[idx].p_cent_distance;
			//integralConst.im = (integralConst1.im - integralConst2.im) * rays[idx].gorden / rays[idx].p_cent_distance;
			integralConst.re = ((integralConst1.re - integralConst2.re) * rays[idx].gorden_re -  (integralConst1.im - integralConst2.im) * rays[idx].gorden_im) / rays[idx].p_cent_distance;
			integralConst.im = ((integralConst1.re - integralConst2.re) * rays[idx].gorden_im -  (integralConst1.im - integralConst2.im) * rays[idx].gorden_re) / rays[idx].p_cent_distance;
			reflect_coeff = config.reflect_coeff_Auto_flag ? d_ReflectCoeff_2(f_ig,rays[idx].Z0) : config.reflect_coeff; //反射系数 姬梓遇20210913
			rays[idx].parameter.p[1] = reflect_coeff * integralConst.re;//反射系数 姬梓遇20210913
			rays[idx].parameter.p[2] = reflect_coeff * integralConst.im;//反射系数 姬梓遇20210913
			
		}
	}
}

//__global__ void statis_sum(float *d_s_sum_re, float *d_sum_re, float *d_s_sum_im, float *d_sum_im, int ig, int j, int maxsize){
	//int idx = threadIdx.x + blockDim.x * blockIdx.x;
	//int tid = threadIdx.x;

	//d_s_sum_re[j * maxsize + ig]= *d_sum_re;
	//d_s_sum_im[j * maxsize + ig]= *d_sum_im;
//}



void scalfuc(RayBeamInfo* rays, int raysBeamNum, int ig, float* d_sum_sre, float* d_sum_sim, ConfigStruct config, Vector* d_center, Axis_slx New_receive_points) {

	dim3 threadSize(16, 1, 1);
	dim3 blockSize(raysBeamNum / 16+ 1, 1, 1);


	transMat << <blockSize, threadSize >> > (rays, raysBeamNum, ig, config, d_center, New_receive_points);//20191218 姬梓遇20210831
	int numItem = 0;



	dim3 THREADSIZE(512, 1, 1);
	dim3 BLOCKSIZE(raysBeamNum / 512 + 1, 1, 1);
	//ReduceAdd.
	reduce_add_sre << <BLOCKSIZE, THREADSIZE >> > (d_sum_sre, rays, raysBeamNum);
	cudaDeviceSynchronize();
	reduce_add_sim << <BLOCKSIZE, THREADSIZE >> > (d_sum_sim, rays , raysBeamNum);
	cudaDeviceSynchronize();
	numItem = BLOCKSIZE.x;
	BLOCKSIZE.x = numItem / THREADSIZE.x + 1;

	while (numItem > 1)
	{
		reduce_sre << <BLOCKSIZE, THREADSIZE >> > (d_sum_sre, d_sum_sre, numItem);
		cudaDeviceSynchronize();
		reduce_sim << <BLOCKSIZE, THREADSIZE >> > (d_sum_sim, d_sum_sim, numItem);
		cudaDeviceSynchronize();
		numItem = BLOCKSIZE.x;
		BLOCKSIZE.x = numItem / THREADSIZE.x + 1;
	}

	//statis_sum << <1 ,1 >> > (d_s_sum_re, d_sum_sre, d_s_sum_im, d_sum_sim, ig - minZeroNum + 200, j, maxsize);
	//float sum_re = 0;
	//float sum_im = 0;
	//HANDLE_ERROR(cudaMemcpy(&sum_re, d_sum_sre, sizeof(float), cudaMemcpyDeviceToHost));
	//HANDLE_ERROR(cudaMemcpy(&sum_im, d_sum_sim, sizeof(float), cudaMemcpyDeviceToHost));

	//comp sum_s = { sum_re, sum_im };

	//return sum_s;



}
