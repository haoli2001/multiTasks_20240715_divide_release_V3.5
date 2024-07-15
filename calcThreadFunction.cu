#include "calcThreadFunction.h"

#include <omp.h>
#include <semaphore.h>
#include <cuda_runtime.h>
#include "handlerror.h"
#include "kd_struct.h"
#include "integral_gpu.h"
#include "raystrace.h"
#include "virtualface_gpu.h"
#include "cudaMallocFree.h"
#include "tree2vector.h"
#include "socketFunctions.h"
#include "simple_time.h"
#include "scalfuc.h"
#include "ReflectCoeff_2.h"
#include "martixMulti.h"
#include <pthread.h>
#define MAX_PIPELINE_CAPICITY 150000000
/**************************
名称：struct PreBlockData
描述：每一个Block中申请的内存指针
***************************/
struct PreBlockData
{
	KD_Node_V *d_tree;
	Prim_Box *d_out_array;
	Element *d_points;
	Triangle *d_triangles;


	Axis direction;
	Direction* d_rays1;
	Direction* d_rays2;
	Square* d_squares1;
	Square* d_squares2;
	RayBeamInfo* d_effrays;
	Vector* d_center;
	Vector* d_axis;
	MatStruct* d_transMat;
	ReimOutput* d_reim;
	float* d_sum_re;
	float* d_sum_im;
	float* d_sum_sre; 
	float* d_sum_sim; 
	int* d_DivRayTubeNum;
	int* d_sum_gmem;
	int* d_sum_Gmem;
	int* d_sum;
	int* d_squares_pred;
	int DivRayTubeNum1st;
	int DivRayTubeNum2nd;
	int DivRayTubeNum3rd;

};

/**************************
名称：struct FreePoint
描述：存储申请的内存指针， 当做参数传递给清理函数cleanup，
      用来在线程被中途杀掉后释放内存
***************************/
struct FreePoint
{
	ConfigStruct config;
	KD_Node_V **h_tree;
	Prim_Box **h_out_array;
	float **pre_triangle_result;
	float **pre_triangle_result_reim;
	PreBlockData **plan;
	int **pre_device_height;
	int **pre_device_width;
	int **height;

	float **e_st_min;
	float **e_fi_max;
	int **divided_num;
    int **divided_width;
    int **divided_height;
    float **divided_st_min;
    float **divided_fi_max;
	ReimOutput*** h_reim;
	comp*** h_TSOfPerTriangle;
	
	DynamicPlane **SubAperturePlane;
	DynamicPlane **AperturePlane;
	BinaryTimeTree ***PreAngelTime;
	
	Element **h_points;	//存储坐标变换后模型point数据的指针
};

/**************************
名称：void cleanup(void* argv)
描述：清理函数，用来在线程被中途杀掉后释放内存
参数：void* argv:   FreePoint 结构体指针
返回值：无
***************************/
void cleanup(void* argv)
{
	FreePoint *fp = (FreePoint*)argv;
	ConfigStruct config = fp->config;
	if((*fp->plan)!=NULL)
	{
		for(int i=0;i<config.card_num;i++)
		{
			HANDLE_ERROR(cudaSetDevice(config.select_device_list[i]));

			FreeOnGPU((*fp->plan)[i].d_rays1, (*fp->plan)[i].d_squares1, (*fp->plan)[i].d_rays2, (*fp->plan)[i].d_squares2, 
				(*fp->plan)[i].d_effrays, (*fp->plan)[i].d_center, (*fp->plan)[i].d_axis, (*fp->plan)[i].d_transMat, (*fp->plan)[i].d_reim, (*fp->plan)[i].d_sum_re, (*fp->plan)[i].d_sum_im, (*fp->plan)[i].d_sum_sre, (*fp->plan)[i].d_sum_sim,
				(*fp->plan)[i].d_DivRayTubeNum, (*fp->plan)[i].d_sum_gmem, (*fp->plan)[i].d_sum_Gmem, (*fp->plan)[i].d_squares_pred);

			free_data((*fp->plan)[i].d_out_array, (*fp->plan)[i].d_tree, (*fp->plan)[i].d_points, (*fp->plan)[i].d_triangles);
		}
	}
	//host Free
	if(*fp->plan!=NULL)
	{
		free(*fp->plan);
	}
	if(*fp->pre_device_height!=NULL)
	{
		free(*fp->pre_device_height);
	}
	if(*fp->pre_device_width!=NULL)
	{
		free(*fp->pre_device_width);
	}
	if(*fp->height!=NULL)
	{
		free(*fp->height);
	}
	if(*fp->e_fi_max!=NULL)
	{
		free(*fp->e_fi_max);
	}
	if(*fp->e_st_min!=NULL)
	{
		free(*fp->e_st_min);
	}
	
	if(*fp->divided_num!=NULL)
	{
		free(*fp->divided_num);
	}
	if(*fp->divided_width!=NULL)
	{
		free(*fp->divided_width);
	}
	if(*fp->divided_height!=NULL)
	{
		free(*fp->divided_height);
	}
	if(*fp->divided_st_min!=NULL)
	{
		free(*fp->divided_st_min);
	}
	if(*fp->divided_fi_max!=NULL)
	{
		free(*fp->divided_fi_max);
	}
	if(*fp->h_out_array!=NULL)
	{
		free(*fp->h_out_array);
	}
	if(*fp->h_tree!=NULL)
	{
		free(*fp->h_tree);
	}
	
	if(*fp->h_reim!=NULL)
	{
		for(int i=0;i<config.card_num;i++)
			free((*fp->h_reim)[i]);
		free(*fp->h_reim);
	}
	if(*fp->h_TSOfPerTriangle!=NULL)
	{
		for(int i=0;i<config.card_num;i++)
			free((*fp->h_TSOfPerTriangle)[i]);
		free(*fp->h_TSOfPerTriangle);
	}
	if(*fp->pre_triangle_result!=NULL)
	{
		free(*fp->pre_triangle_result);
	}
	if(*fp->pre_triangle_result_reim!=NULL)
	{
		free(*fp->pre_triangle_result_reim);
	}
	if(*fp->SubAperturePlane!=NULL)
	{
		free(*fp->SubAperturePlane);
	}
	if(*fp->AperturePlane!=NULL)
	{
		free(*fp->AperturePlane);
	}
	if(*fp->PreAngelTime!=NULL)
	{
		for(int i=0;i < 2 * config.card_num - 1;i++)
			free((*fp->PreAngelTime)[i]);
		free(*fp->PreAngelTime);
	}
	if(*fp->h_points!=NULL)
	{
		free(*fp->h_points);
	}
}

/**************************
名称：void* calcThreadFunction(void *argv)
描述：计算线程函数
参数：void* argv:   传递给计算线程的参数
返回值：无
***************************/
void* calcThreadFunction(void *argv)
{
    //计算线程参数
	CalcInfo calcInfo = *(CalcInfo*)argv;
	
	int socketClient = calcInfo.socket;     
	KD_Node_V *h_tree = NULL;
	Prim_Box *h_out_array = NULL;
	Element *h_points = NULL;
	Element *h_points_old = NULL;
	Triangle *h_triangles = NULL;
	float *pre_triangle_result = NULL;
	float *pre_triangle_result_reim = NULL;
	PreBlockData *plan = NULL;
	int *pre_device_height = NULL;
	int *pre_device_width = NULL;
	int *height = NULL;
	int *divided_num = NULL;
	int *divided_width = NULL;
	int *divided_height = NULL;
	float *divided_st_min = NULL;
	float *divided_fi_max = NULL;
	float *e_st_min = NULL;
	float *e_fi_max = NULL;
	ReimOutput** h_reim = NULL;
	comp** h_TSOfPerTriangle = NULL;
	DynamicPlane *SubAperturePlane = NULL;
	DynamicPlane *AperturePlane = NULL;
	BinaryTimeTree **PreAngelTime = NULL;
	
    //中途退出需要释放的内存结构体
	FreePoint freePt;
	freePt.config= calcInfo.config;
	freePt.h_tree=&h_tree;
	freePt.h_out_array=&h_out_array;
	freePt.pre_triangle_result=&pre_triangle_result;
	freePt.pre_triangle_result_reim=&pre_triangle_result_reim;
	freePt.plan=&plan;
	freePt.pre_device_height=&pre_device_height;
	freePt.pre_device_width=&pre_device_width;
	freePt.height=&height;

	freePt.e_st_min=&e_st_min;
	freePt.e_fi_max=&e_fi_max;
	freePt.divided_num=&divided_num;
	freePt.divided_width=&divided_width;
	freePt.divided_height=&divided_height;
	freePt.divided_st_min=&divided_st_min;
	freePt.divided_fi_max=&divided_fi_max;
	freePt.h_reim=&h_reim;
	freePt.h_TSOfPerTriangle=&h_TSOfPerTriangle;
	
	freePt.SubAperturePlane=&SubAperturePlane;
	freePt.AperturePlane=&AperturePlane;
	freePt.PreAngelTime=&PreAngelTime;
	freePt.h_points = &h_points;

    //设置退出线程时调用cleanup清理内存
	pthread_cleanup_push(cleanup,&freePt);
    //允许退出线程 
	pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
    //下一个监视点取消    
	pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED, NULL);  

	int h_tree_length;
	int h_out_array_length;
	int h_points_length;
	int h_triangle_length;
	
	//每一个三角形积分结果，用于伪彩图显示
	pre_triangle_result = (float*)malloc(sizeof(float)*calcInfo.triangles_length);
	memset(pre_triangle_result,0,sizeof(float)*calcInfo.triangles_length);
	pre_triangle_result_reim = (float*)malloc(2*sizeof(float)*calcInfo.triangles_length);
	memset(pre_triangle_result_reim,0,2*sizeof(float)*calcInfo.triangles_length);
    
	//根据船体晃动的三个角度，生成坐标旋转矩阵，旋转模型坐标
	//并根据用户的输入水线深度，计算模型坐标中的绝对水线z坐标
	float rolX[9] = {1,0,0,0,cos(calcInfo.config.sway_phi),sin(calcInfo.config.sway_phi),0,-sin(calcInfo.config.sway_phi),cos(calcInfo.config.sway_phi)};
	float rolY[9] = {cos(calcInfo.config.sway_theta),0,-sin(calcInfo.config.sway_theta),0,1,0,sin(calcInfo.config.sway_theta),0,cos(calcInfo.config.sway_theta)};
	float rolZ[9] = {cos(calcInfo.config.sway_psi),sin(calcInfo.config.sway_psi),0,-sin(calcInfo.config.sway_psi),cos(calcInfo.config.sway_psi),0,0,0,1};
	float Change_Corrd_Matrix[9];
	float tempMatrix[9];
	martixMulti(rolX, 3, 3, rolZ, 3, 3, tempMatrix);
	martixMulti(tempMatrix, 3, 3, rolY, 3, 3, Change_Corrd_Matrix);
	float abs_waterLine_axis = INT_MAX;
	h_points_old = calcInfo.points;
	h_points_length = calcInfo.points_length;
	h_points = (Element*)malloc(sizeof(Element) * h_points_length);
	for(int i=0;i<h_points_length;i++)
	{
		if(h_points_old[i].point[2] < abs_waterLine_axis)
			abs_waterLine_axis = h_points_old[i].point[2];
		martixMulti(Change_Corrd_Matrix, 3, 3, h_points_old[i].point, 3, 1, h_points[i].point);
	}
	abs_waterLine_axis = abs_waterLine_axis + calcInfo.config.water_line;
	printf("\n模型坐标转换完毕，水线深度：%f，模型中水线z坐标：%f\n",calcInfo.config.water_line,abs_waterLine_axis);


	//记录KDTree 建树时间
	simple_time simpleTimeKDTree;
	float timeKDTree;

    //建立KDTree
	printf("\n\n*************KD_Tree build begin*************\n");
	simpleTimeKDTree.Time_Start();
	{
		int Out_Arr_Length = 0;
		int index = 0;
		
		KD_Node *Root_Node = (KD_Node*)malloc(sizeof(KD_Node));
		Prim_Box *Triangles_Box = (Prim_Box*)malloc(calcInfo.triangles_length * sizeof(Prim_Box));

		Triangle *Triangles;
		Element *Points;
		Prim_Box *Out_Array;
		KD_Infomation KD_Info;

		Out_Array= (Prim_Box *)malloc( 4 * calcInfo.triangles_length * sizeof(Prim_Box));

		Sort_Box(calcInfo.points, calcInfo.triangles, Triangles_Box, calcInfo.triangles_length);
		memset(Root_Node, 0, sizeof(struct KD_Node));//��ʼ�����ڵ�Ϊ0
		memset(&KD_Info, 0, sizeof(struct KD_Infomation));//��ʼ��ͳ����ϢΪ0
		KD_Node_init(Root_Node, Triangles_Box, calcInfo.triangles_length);//��������Ķ�����Ϣ�����г�ʼ����ֵ


		Build_BigNode(Root_Node, Triangles_Box, Out_Array, calcInfo.triangles_length, &Out_Arr_Length);
		Optimization_Rope(Root_Node);//����ֵ���Ż�

        //将KDTree之前的树状存储结构转换为连续存储的数组形式
		KD_Node_V *root_v = tree2vector(Root_Node, &index);
		
        //销毁树状结构的KDTree
		Destroy_Tree(Root_Node);
		
		free(Triangles_Box);
	
		h_out_array = Out_Array;
		h_out_array_length = Out_Arr_Length;
		h_tree = root_v;
		h_tree_length = index;
		h_triangles = calcInfo.triangles;
		h_triangle_length = calcInfo.triangles_length;
	}
	timeKDTree = simpleTimeKDTree.Time_End();

	//有时KDTree的length计算会出错，使得该值非常大，因而增加此判断
	if(h_tree_length > 100000000)
		exit(0);

	printf("发送KDTree,建树时间，节点个数等信息\n");
	{
		struct KDTreeInfo
		{
			float time;
			int length;
		};
		KDTreeInfo kdtreeinfo;
		kdtreeinfo.time = timeKDTree;
		kdtreeinfo.length = h_tree_length;
		Frame frame;
		strcpy(frame.command, "KDTreeTime");
		frame.length = sizeof(KDTreeInfo);
		memcpy(frame.data, &kdtreeinfo, sizeof(KDTreeInfo));
		send_frame(socketClient, (char*)&frame, sizeof(Frame));
	}
    
    //线程退出的监视点
	pthread_testcancel(); 
	pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, NULL);
    
	printf("发送KDTree数据！\n");
	{
		int sendedLength = 0;
		Frame frame;

		while (true)
		{
			strcpy(frame.command, "KDTreeDate");
			if (sendedLength + 1024 < h_tree_length*sizeof(KD_Node_V))
			{
				memcpy(frame.data, (char*)h_tree + sendedLength, 1024);
				frame.length = 1024;
				send_frame(socketClient, (char*)&frame, sizeof(Frame));
				sendedLength += 1024;
			}
			else
			{
				memcpy(frame.data, (char*)h_tree + sendedLength, h_tree_length*sizeof(KD_Node_V) - sendedLength);
				frame.length = h_tree_length*sizeof(KD_Node_V) - sendedLength;
				send_frame(socketClient, (char*)&frame, sizeof(Frame));
				break;
			}
		}
	}
	
	pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
    //线程退出的监视点
	pthread_testcancel(); 
	pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, NULL);

	plan = (PreBlockData*)malloc(sizeof(PreBlockData)*calcInfo.config.card_num);
	memset(plan, 0, sizeof(PreBlockData)*calcInfo.config.card_num);

	//int sumAngleNum = calcInfo.config.end_alpha - calcInfo.config.start_alpha + 1;
	//保证初始角度不为0时仍能正确计算,  220304 jzy
	int sumAngleNum = calcInfo.config.end_alpha + 1;
	//每一个卡上的宽
	pre_device_height = (int*)malloc(sizeof(int)*sumAngleNum);
	//每一个卡上的高
	pre_device_width = (int*)malloc(sizeof(int)*sumAngleNum);
	divided_num = (int*)malloc(sizeof(int)*sumAngleNum);
	height = (int*)malloc(sizeof(int)*sumAngleNum);

	e_st_min = (float*)malloc(sizeof(float)*sumAngleNum);
	e_fi_max = (float*)malloc(sizeof(float)*sumAngleNum);
	int d_width_max = -1;
	int d_height_max = -1;
    float far_dis = calcInfo.config.far_distance;


	//多卡上的malloc 和树，顶点，三角面片，包围盒,使用OpenMp实现多线程
	omp_set_num_threads(calcInfo.config.card_num);
	#pragma omp parallel
	{
		int i = omp_get_thread_num();
		cudaSetDevice(calcInfo.config.select_device_list[i]);

		HANDLE_ERROR(cudaMalloc((void**)&plan[i].d_tree, sizeof(KD_Node_V)*h_tree_length));
		HANDLE_ERROR(cudaMalloc((void**)&plan[i].d_out_array, sizeof(Prim_Box)*h_out_array_length));
		HANDLE_ERROR(cudaMalloc((void**)&plan[i].d_points, sizeof(Element)*h_points_length));
		HANDLE_ERROR(cudaMalloc((void**)&plan[i].d_triangles, sizeof(Triangle)*h_triangle_length));
		cudaMemcpy(plan[i].d_tree, h_tree, h_tree_length*sizeof(KD_Node_V), cudaMemcpyHostToDevice);
		cudaMemcpy(plan[i].d_out_array, h_out_array, h_out_array_length*sizeof(Prim_Box), cudaMemcpyHostToDevice);
		cudaMemcpy(plan[i].d_points, h_points, h_points_length*sizeof(Element), cudaMemcpyHostToDevice);
		cudaMemcpy(plan[i].d_triangles, h_triangles, h_triangle_length*sizeof(Triangle), cudaMemcpyHostToDevice);
	}

	float fai_angle = calcInfo.config.start_beta;
	//获取每一度虚拟孔径面所需宽高，通过先得到每一度的宽高，按最大的开辟空间，减少了每一度扫描的空间开辟时间
	for (int i = calcInfo.config.start_alpha; i <= calcInfo.config.end_alpha; i++)
	{
		//注意：此处angle没有考虑浮点角度！另：pre_device_height无用	jzy 2023.4.4
		float angle = i;
		getWidthHeight(far_dis,h_points, h_points_length, &pre_device_width[i], &height[i], si_angle, angle, calcInfo.config.pipe_size*calcInfo.config.wave_length, &pre_device_height[i], &e_st_min[i], &e_fi_max[i], calcInfo.config.card_num, &divided_num[i], MAX_PIPELINE_CAPICITY);
	}
	//printf("d_height_max=%d, d_width_max=%d\n", d_height_max, d_width_max);
	d_height_max = 1975;
	d_width_max = 10000;

	pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
    //线程退出的监视点	
	pthread_testcancel(); 
    pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, NULL);

	//根据最大的宽高开辟空间
	omp_set_dynamic(0);
	omp_set_num_threads(calcInfo.config.card_num);  // create as many CPU threads as there are CUDA devices
	#pragma omp parallel
	{

		int i = omp_get_thread_num();
		cudaSetDevice(calcInfo.config.select_device_list[i]);

		MallocOnGPU(d_width_max, d_height_max, &plan[i].d_rays1, &plan[i].d_squares1, &plan[i].d_rays2, &plan[i].d_squares2,
			&plan[i].d_effrays, &plan[i].d_center, &plan[i].d_axis, &plan[i].d_transMat, &plan[i].d_reim, &plan[i].d_sum_re, &plan[i].d_sum_im, &plan[i].d_sum_sre, &plan[i].d_sum_sim,
			&plan[i].d_DivRayTubeNum, &plan[i].d_sum_gmem, &plan[i].d_sum_Gmem, &plan[i].d_squares_pred);
	}

	h_reim = (ReimOutput**)malloc(sizeof(ReimOutput*)*calcInfo.config.card_num);
	for (int t = 0; t < calcInfo.config.card_num; t++)
	{
		h_reim[t] = (ReimOutput*)malloc(sizeof(ReimOutput) * d_width_max * d_height_max);
		memset(h_reim[t], 0, sizeof(ReimOutput) * d_width_max * d_height_max);
	}

	h_TSOfPerTriangle = (comp**)malloc(sizeof(comp*) * calcInfo.config.card_num);
	for (int p = 0; p < calcInfo.config.card_num; p++)
	{
		h_TSOfPerTriangle[p] = (comp*)malloc(sizeof(comp) * calcInfo.triangles_length);
		memset(h_TSOfPerTriangle[p], 0, sizeof(comp) * calcInfo.triangles_length);
	}
	
	int NodeNum = 2 * calcInfo.config.card_num - 1;
	SubAperturePlane = (DynamicPlane*)malloc(calcInfo.config.card_num * sizeof(DynamicPlane));
	AperturePlane = (DynamicPlane*)malloc(NodeNum * sizeof(DynamicPlane));
	PreAngelTime = (BinaryTimeTree**)malloc(NodeNum * sizeof(BinaryTimeTree*));
	for (int i = 0; i < NodeNum; i++)
	{
		PreAngelTime[i] = (BinaryTimeTree*)malloc(sizeof(BinaryTimeTree));
		memset(PreAngelTime[i], 0, sizeof(BinaryTimeTree));
	}

	float CSpeed = 1500.0;//20200919
	float pi = 3.1415926535898;
	float D2R = pi / 180.0;
	float fend = calcInfo.config.time_end_frequency;//姬梓遇20210831
	float fbeg = calcInfo.config.time_start_frequency;////姬梓遇20210831
	float Tao = calcInfo.config.tao/1000;  //原单位为ms，此处转为s 姬梓遇20210831
	float fs = calcInfo.config.sampling_frequency;//caiyanglv 姬梓遇20210831
	float velocity1 = calcInfo.config.velocity1;//鱼雷速度1 姬梓遇
	float velocity2 = calcInfo.config.velocity2;//目标速度1 姬梓遇
    float velocity12 = velocity1 - velocity2;//相对投影速度 LV 20220720
    //float velocity12 = velocity1 * cos(calcInfo.config.H_axis_angle * D2R)  -  velocity2 * cos(calcInfo.config.? + calcInfo.config.H_axis_angle * D2R);//相对投影速度 LV 20220720
	float band = fend - fbeg;
	float KO = band / Tao;//wangying	
    float dopGene = (CSpeed + velocity12) / (CSpeed - velocity12);//多普勒压缩因子
	if(fabs(band) > 0.001)  //宽带信号则考虑脉冲宽度压缩
		Tao = Tao / dopGene;
	fbeg = fbeg * dopGene;//调制后起始频率
	fend = fbeg + band * dopGene * dopGene;   // LV20221229 原有计算，错误： fend = fbeg + band;

	int taosize = int(fs * Tao);
	int maxsize;
	int totalsize= int(calcInfo.config.sampling_width * fs);//20210402

	//FILE* fileresult_1200 = fopen("result_1200", "w");//wangying
	FILE* fileresult_1200;
    FILE* fileresult_TS = fopen("TS.txt", "w");

	char nameout[50];
	float* result_1200_re;
	//memset(result_1200_re, 0, maxsize * sizeof(float));
	float* result_1200_im;
	//memset(result_1200_im, 0, maxsize * sizeof(float));//20200920
	int minZeroNum=INT_MAX,maxZeroNum=INT_MIN;//20210308

	float* result_1200;//20210402
	comp result={0,0};
	float* s_sum_re;//20200919 面积积分结果
	float* s_sum_im;//20200919 面积积分结果

	float* d_m_sum_re[calcInfo.config.card_num];//20210331
	float* d_m_sum_im[calcInfo.config.card_num];//20210331


	//生成多点接收阵列20220527jzy
	int sumnum = calcInfo.config.recvPointsNum;
	float d = 0.026;
	float theta = 0;
	float psi = 0;
	float fai = 0;
	//Axis_slx* receive_points = (Axis_slx*)malloc(sizeof(Axis_slx) * sumnum);
	Axis_slx* receive_points = calcInfo.recvPoints;
	Axis_slx* New_receive_points = (Axis_slx*)malloc(sizeof(Axis_slx) * sumnum);
	printf("雷体坐标系下各坐标\n");
	for(int idx=0;idx<sumnum;idx++)
	{		
            printf("%f %f %f\n",receive_points[idx].p[0],receive_points[idx].p[1],receive_points[idx].p[2]);
	}
	//生成多点接收阵列结束20220527jzy




	//循环计算每一个频率和角度
	float start_alpha;
	if(fabs(calcInfo.config.continue_alpha + 1) < 0.00001)
		start_alpha = calcInfo.config.start_alpha;
	else
		start_alpha = calcInfo.config.continue_alpha;
	

	for(int recv_index = calcInfo.config.recvPointsStartIdx; recv_index < sumnum; recv_index ++)
	{
		float angle = start_alpha;
		for (float f = calcInfo.config.start_frequency; f <= calcInfo.config.end_frequency; f += 0.1)
		{
			for (int i = start_alpha; i <= calcInfo.config.end_alpha; i++)
			{
				float lmd = 1.5 / f;          // 声波波长
				float time_all=0;
				result.re=0;
				result.im=0;
				memset(result_1200, 0, totalsize * sizeof(float));

				divided_width = (int*)malloc(sizeof(int) * divided_num[i]);
				divided_height = (int*)malloc(sizeof(int) * divided_num[i]);
				divided_st_min = (float*)malloc(sizeof(float) * divided_num[i]);
				divided_fi_max = (float*)malloc(sizeof(float) * divided_num[i]);
				//分割虚拟孔径面
				divide_module_virtualface(divided_num[i], calcInfo.config.pipe_size*calcInfo.config.wave_length, sumAngleNum, pre_device_width[i], height[i], e_st_min[i], e_fi_max[i], divided_width, divided_height, divided_st_min, divided_fi_max);






            	//多卡的时间统计 同时启动
            	simple_time *runSimpleTime = new simple_time[calcInfo.config.card_num];
           		float *calcTime = (float*)malloc(sizeof(float)*calcInfo.config.card_num);
           		memset(calcTime,0,sizeof(float)*calcInfo.config.card_num);
				// = fopen("i", "w");//wangying
				sprintf(nameout,"./data/%.2f_recv%d.txt",angle,recv_index);
				fileresult_1200=fopen(nameout,"w");//20200924
		
				pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
            	//线程退出的监视点
				pthread_testcancel(); 
				//避免线程在其他地方cancel，以保证上位机能正常暂停计算 2022.3.24 jzy
				pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, NULL);
            	//先将数据清零，避免上一度的结果影响
				for (int p = 0; p < calcInfo.config.card_num; p++)
				{
					memset(h_TSOfPerTriangle[p], 0, sizeof(comp) * calcInfo.triangles_length);
				}
				comp *rnt_sum = (comp *)malloc(sizeof(comp)*calcInfo.config.card_num);
				memset(rnt_sum, 0,sizeof(comp)*calcInfo.config.card_num);
				float lmd = 1.5 / f;          // 声波波长


				//多点接收阵列坐标转换到目标坐标系
				float R0 = far_dis;
				float alpha = angle * D2R;
				float beta = fai_angle * D2R;
				//float P0_B_A = R0 * sin(beta) * cos(alpha)
				float H_Axis_Angle = calcInfo.config.H_axis_angle * D2R;
				float V_Axis_Angle = calcInfo.config.V_axis_angle * D2R;
		
				printf("H_Axis_Angle = %f\n",H_Axis_Angle);	//调试用
				if((alpha >=0* D2R)&&(alpha <90*D2R))
				{
					H_Axis_Angle = H_Axis_Angle;
				}
				else if((alpha >=90* D2R)&&(alpha <=180*D2R))
				{
					H_Axis_Angle = -H_Axis_Angle;
				}
				else
				{}

				if((beta  >=0* D2R)&&(beta  <90*D2R))
				{
					V_Axis_Angle = -V_Axis_Angle;
				}
				else if((beta  >90* D2R)&&(beta  <=180*D2R))
				{
					V_Axis_Angle = V_Axis_Angle;
				}
				else
				{}

				psi = (beta + V_Axis_Angle - pi/2);
				theta = -(pi - (alpha - H_Axis_Angle));
				//生成坐标转换矩阵


				//float rolX[9] = {1,0,0,0,cos(fai),sin(fai),0,-sin(fai),cos(fai)};
				//float rolY[9] = {cos(theta),0,-sin(theta),0,1,0,sin(theta),0,cos(theta)};
				//float rolZ[9] = {cos(psi),sin(psi),0,-sin(psi),cos(psi),0,0,0,1};


				float rolX[9] = {1,0,0,0,cos(fai),sin(fai),0,-sin(fai),cos(fai)};
				float rolY[9] = {cos(theta),0,-sin(theta),0,1,0,sin(theta),0,cos(theta)};
				float rolZ[9] = {cos(psi),sin(psi),0,-sin(psi),cos(psi),0,0,0,1};

				float Change_Corrd_Matrix[9];
				float tempMatrix[9];
				martixMulti(rolX, 3, 3, rolZ, 3, 3, tempMatrix);
				martixMulti(tempMatrix, 3, 3, rolY, 3, 3, Change_Corrd_Matrix);
				int jj;
				/*
				printf("fai:%f\n",fai);
				printf("theta:%f\n",theta);
				printf("psi:%f\n",psi);

			
				printf("rolX旋转矩阵各坐标\n");
				for(jj = 0;jj<9;jj++)
				{
					printf("%f\n",rolX[jj]);
				}
				printf("rolY旋转矩阵各坐标\n");
				for(jj = 0;jj<9;jj++)
				{
					printf("%f\n",rolY[jj]);
				}
				for(jj = 0;jj<9;jj++)
				{
					printf("%f\n",rolZ[jj]);
				}
				printf("旋转矩阵各坐标\n");
				for(jj = 0;jj<9;jj++)
				{
					printf("%f\n",Change_Corrd_Matrix[jj]);
				}
				printf("目标坐标系下基阵几何中心坐标\n");	
				printf("%f %f %f\n",R0 * sin(beta) * cos(alpha),R0 * cos(beta),R0 * sin(beta) * sin(alpha));
				*/
				printf("目标坐标系下各坐标\n");	
				for(int idx=0;idx<sumnum;idx++)
				{
					martixMulti(Change_Corrd_Matrix, 3, 3, receive_points[idx].p, 3, 1, New_receive_points[idx].p);
					printf("%f %f %f\n",New_receive_points[idx].p[0],New_receive_points[idx].p[1],New_receive_points[idx].p[2]);
					New_receive_points[idx].p[0] += R0 * sin(beta) * cos(alpha);
					New_receive_points[idx].p[1] += R0 * cos(beta);
					New_receive_points[idx].p[2] += R0 * sin(beta) * sin(alpha);

					//互换y、z坐标，保证一致
					float temp_index = 0;
					temp_index = New_receive_points[idx].p[1];
					New_receive_points[idx].p[1] = New_receive_points[idx].p[2];
					New_receive_points[idx].p[2] = temp_index;

				}
			
				//多点接收阵列坐标转换到目标坐标系结束


				//	float *s_sum_re;//20200919 面积积分结果
				//float* s_sum_im;//20200919 面积积分结果

				// 动态生成子孔径面边界信息
				ConstructVirtualFace(AperturePlane, SubAperturePlane, PreAngelTime, calcInfo.config.card_num, e_st_min[i], e_fi_max[i], pre_device_width[i], height[i], calcInfo.config.pipe_size * calcInfo.config.wave_length);
			
				// 开始计时
				for(int index=0;index<calcInfo.config.card_num;index++)
            	{
               	 	runSimpleTime[index].Time_Start();
            	}
			
				//多卡并行计算
            	simple_time simple_time_TS_compute;
            	simple_time_TS_compute.Time_Start();
				omp_set_num_threads(calcInfo.config.card_num);  // create as many CPU threads as there are CUDA device
				#pragma omp parallel
				{
				int j = omp_get_thread_num();                        //目前线程id,即卡id
				int num_threads = omp_get_num_threads();             //获取卡的数量
				cudaSetDevice(calcInfo.config.select_device_list[j]);
				//cudaDeviceEnablePeerAccess(j, 0);

				Radius direction_radius = { far_dis, fai_angle, angle }; //snw 接入远场距离
				//printf("direction_radius={%f,%f,%f}\n",direction_radius.Xr,direction_radius.Yst,direction_radius.Zfi);
				plan[j].direction = dSphericaltoRectangular(direction_radius);//等相位面的法向量

				MemsetOnGPU(d_width_max, d_height_max, &plan[j].d_rays1, &plan[j].d_squares1, &plan[j].d_rays2, &plan[j].d_squares2,
					&plan[j].d_effrays, &plan[j].d_center, &plan[j].d_axis, &plan[j].d_transMat, &plan[j].d_reim, &plan[j].d_sum_re, &plan[j].d_sum_im, &plan[j].d_sum_sre, &plan[j].d_sum_sim,
					&plan[j].d_DivRayTubeNum, &plan[j].d_sum_gmem, &plan[j].d_sum_Gmem, &plan[j].d_squares_pred);
					
				//printf("GPU%d: width=%d, height=%d, st_min=%e, fi_max=%e\n", j, SubAperturePlane[j].width, SubAperturePlane[j].height, SubAperturePlane[j].st_min, SubAperturePlane[j].fi_max);

                //创建虚拟孔径面
				create_virtualface_gpu(plan[j].d_rays1, plan[j].d_squares1, SubAperturePlane[j].width, SubAperturePlane[j].height,
					calcInfo.config.pipe_size * calcInfo.config.wave_length, direction_radius, SubAperturePlane[j].st_min, SubAperturePlane[j].fi_max, j);
                
                //射线追踪
				allraystrace_v2(plan[j].d_rays1, plan[j].d_squares1, SubAperturePlane[j].width, SubAperturePlane[j].height,
					plan[j].d_tree, plan[j].d_out_array, plan[j].d_points, plan[j].d_triangles, plan[j].d_DivRayTubeNum,
					&(plan[j].DivRayTubeNum1st), plan[j].d_sum_gmem, plan[j].d_sum_Gmem, plan[j].d_squares_pred, plan[j].direction, angle, abs_waterLine_axis);
                
                //声场积分
				RayBeamInfo* c_effrays = (RayBeamInfo*)malloc(SubAperturePlane[j].width * SubAperturePlane[j].height  * sizeof(RayBeamInfo));

				comp sum = sound_field_integral_gpu(plan[j].d_rays1, plan[j].d_squares1, lmd, SubAperturePlane[j].width * SubAperturePlane[j].height, plan[j].d_effrays, plan[j].d_center,
					plan[j].d_axis, plan[j].d_transMat, plan[j].d_reim, plan[j].d_sum_re, plan[j].d_sum_im, fai_angle, angle, &(calcInfo.config));//20210831姬梓遇
                calcTime[j]+=runSimpleTime[j].Time_End();

				HANDLE_ERROR(cudaMemcpy(c_effrays, plan[j].d_effrays, SubAperturePlane[j].width * SubAperturePlane[j].height * sizeof(RayBeamInfo), cudaMemcpyDeviceToHost)); 


				//copy出各声线的d_center,用于下面计算时域积分的maxsize	20220610
				Vector* c_center = (Vector*)malloc(SubAperturePlane[j].width * SubAperturePlane[j].height * sizeof(Vector));
				HANDLE_ERROR(cudaMemcpy(c_center, plan[j].d_center, SubAperturePlane[j].width * SubAperturePlane[j].height * sizeof(Vector), cudaMemcpyDeviceToHost)); 

				for(int idx=0;idx<SubAperturePlane[j].width * SubAperturePlane[j].height;idx++)
				{
					float recv_p_cent_distance = sqrt(pow(c_center[idx].x - New_receive_points[recv_index].p[0], 2) + pow(c_center[idx].y - New_receive_points[recv_index].p[1], 2) + pow(c_center[idx].z - New_receive_points[recv_index].p[2], 2));
					int pendZeroNum = int( (recv_p_cent_distance + c_effrays[idx].p_cent_distance) / CSpeed * fs) ;//wangying  snw:pendZeroNum从2倍距离算起
					if(pendZeroNum!=0&&pendZeroNum<minZeroNum) minZeroNum=pendZeroNum;
					if(pendZeroNum>maxZeroNum) maxZeroNum=pendZeroNum;
				}//20210308
				free(c_effrays);
        		free(c_center);	   
                //拷贝出积分结果，用来伪彩图显示
				cudaMemcpy(h_reim[j], plan[j].d_reim, d_width_max * d_height_max * sizeof(ReimOutput), cudaMemcpyDeviceToHost);
				for (int n = 0; n < d_width_max * d_height_max; n++)
				{
					if (h_reim[j][n].triangle_index >= 0)
					{
						h_TSOfPerTriangle[j][h_reim[j][n].triangle_index].re += h_reim[j][n].re;
						h_TSOfPerTriangle[j][h_reim[j][n].triangle_index].im += h_reim[j][n].im;
					}
				}
				float reflect_coeff = calcInfo.config.reflect_coeff_Auto_flag ? ReflectCoeff_2(f,i) : calcInfo.config.reflect_coeff; //计算反射系数 姬梓遇
				sum.re = sum.re * reflect_coeff;//积分结果乘反射系数
				sum.im = sum.im * reflect_coeff;//积分结果乘反射系数
				rnt_sum[j].im = sum.im;
				rnt_sum[j].re = sum.re;

				maxsize = maxZeroNum + taosize - minZeroNum + 400;

				printf("maxsize, minZeroNum, maxZeroNum+taosize: %d, %d, %d\n",maxsize, minZeroNum, maxZeroNum+taosize);

				//float *d_m_sum_re,*d_m_sum_im;
				
				//20210331
				HANDLE_ERROR(cudaMalloc((void** )&d_m_sum_re[j], maxsize * sizeof(float)));
				HANDLE_ERROR(cudaMalloc((void** )&d_m_sum_im[j], maxsize * sizeof(float)));
				
				for (int ig = 0; ig <  maxsize; ig++) 
				{//wangying 20210308
					float tao0 = ig * (1 / fs);
					float wavelength = 2*PI*(fbeg+KO*tao0)/CSpeed;
					scalfuc(plan[j].d_effrays, SubAperturePlane[j].width * SubAperturePlane[j].height, ig + minZeroNum - 200, d_m_sum_re[j] + ig, d_m_sum_im[j] + ig, calcInfo.config, plan[j].d_center, New_receive_points[recv_index]);//20210331 20210831姬梓遇
                    //scalfuc(plan[j].d_effrays, SubAperturePlane[j].width * SubAperturePlane[j].height, ig + minZeroNum - 200, d_m_sum_re[j] + ig, d_m_sum_im[j] + ig, calcInfo.config, d_one_beam_result[j] + ig);//20210904snw   
					//printf("ig * sizeof(float): %d\n", ig * sizeof(float));

					//plan[j].d_m_sum_re[ig]=*plan[j].d_sum_sre;
					//plan[j].d_m_sum_im[ig]=*plan[j].d_sum_sim;
				}

				

			}
		
			s_sum_re = (float* )malloc(maxsize * sizeof(float) * calcInfo.config.card_num);//20200919 面积积分结果
			s_sum_im = (float* )malloc(maxsize * sizeof(float) * calcInfo.config.card_num);//20200919 面积积分结果	
			//20210331
			omp_set_num_threads(calcInfo.config.card_num);  // create as many CPU threads as there are CUDA devices
			#pragma omp parallel
			{
				int j = omp_get_thread_num();
				HANDLE_ERROR(cudaSetDevice(calcInfo.config.select_device_list[j]));
		
				HANDLE_ERROR(cudaMemcpy(s_sum_re+ j * maxsize, d_m_sum_re[j], maxsize * sizeof(float), cudaMemcpyDeviceToHost));
				HANDLE_ERROR(cudaMemcpy(s_sum_im+ j * maxsize, d_m_sum_im[j], maxsize * sizeof(float), cudaMemcpyDeviceToHost));
				cudaFree(d_m_sum_re[j]);
				cudaFree(d_m_sum_im[j]);
			}	
			
            //simple_time simple_time_TS_compute;
            //simple_time_TS_compute.Time_Start();
            //合并每一块卡的积分结果
			comp all_rnt_sum = {0,0};
			for(int index=0;index<calcInfo.config.card_num;index++)
			{
				all_rnt_sum.re += rnt_sum[index].re;
				all_rnt_sum.im += rnt_sum[index].im;
			}

			//memset(result_1200_re, 0, maxsize * sizeof(float));
			//memset(result_1200_im, 0, maxsize * sizeof(float));
			
			free(rnt_sum);

		 	result_1200_re = (float*)malloc(maxsize * sizeof(float));
			memset(result_1200_re, 0, maxsize * sizeof(float));
			result_1200_im = (float*)malloc(maxsize * sizeof(float));
			memset(result_1200_im, 0, maxsize * sizeof(float));//20200920

			for (int ig =  0; ig <  maxsize; ig++)//wangying 20210308 
			{
				float all_s_sum_re = 0;
				float all_s_sum_im = 0;
				for (int index = 0; index < calcInfo.config.card_num; index++)
				{
					all_s_sum_re += s_sum_re[index * maxsize + ig ];
					all_s_sum_im += s_sum_im[index * maxsize + ig ];
				}
				result_1200_re[ig] = all_s_sum_re;
				result_1200_im[ig] = all_s_sum_im;
			}



			free(s_sum_re);
			free(s_sum_im);


			result_1200 = (float*)malloc( totalsize* sizeof(float));//20210402
			memset(result_1200, 0, totalsize * sizeof(float));
			for(int ig =  0; ig <  maxsize; ig++){
				result_1200[ig+minZeroNum-200]=result_1200_re[ig];//20220719 考虑-200
				}




			for (int k = 0; k < totalsize; k++)//20210402
			{
				//fprintf(fileresult_1200, "%d	%d	%e	%e\n", i,k, result_1200_re[k], result_1200_im[k]);
				fprintf(fileresult_1200, "%e\n", result_1200[k]);
			}
			
			fclose(fileresult_1200);
            //计算TS值
			float result = TS_compute(far_dis,all_rnt_sum, lmd);
            
            float time_TS_compute = simple_time_TS_compute.Time_End();
			
			float maxcalcTime = -1;
            for(int index=0;index<calcInfo.config.card_num;index++)
			{
				if(calcTime[index]>maxcalcTime)
                    maxcalcTime = calcTime[index];
			}
            
            float time_all = maxcalcTime+time_TS_compute;
            
            
            //计算每一个三角面片的积分结果，用于伪彩图显示
			memset(pre_triangle_result,0,sizeof(float)*calcInfo.triangles_length);
			for (int j = 0; j < calcInfo.config.card_num; j++)
			{
				for (int index = 0; index < calcInfo.triangles_length; index++)
				{
					pre_triangle_result_reim[2*index] += h_TSOfPerTriangle[j][index].re;
					pre_triangle_result_reim[2*index+1] += h_TSOfPerTriangle[j][index].im;
					if (h_TSOfPerTriangle[j][index].re != 0 || h_TSOfPerTriangle[j][index].im != 0)
					{				
						pre_triangle_result[index] += sqrt(h_TSOfPerTriangle[j][index].re*h_TSOfPerTriangle[j][index].re + h_TSOfPerTriangle[j][index].im*h_TSOfPerTriangle[j][index].im);						
					}
				}
			}
			
			// 将各卡计算时间存储成二叉树结构，用于调整下一度子孔径面划分
			PreAngelTime = ConstructTimeTree(AperturePlane, calcTime, calcInfo.config.card_num);
		
			CalcResult calcResult;
			calcResult.angle = angle;
			calcResult.calc_time = time_all;
			calcResult.freq = f;
			calcResult.raysnum = pre_device_width[i] * height[i] + (pre_device_width[i] + 1)*(height[i] + 1);
			calcResult.squarenum = pre_device_width[i] * height[i];
			calcResult.TS = result;
			calcResult.height = height[i];
			calcResult.width = pre_device_width[i];
			calcResult.recvIdx = recv_index;

            fprintf(fileresult_TS,"%e\n",calcResult.TS);  //snw

			Frame frame;
			strcpy(frame.command, "CalcResult");
			frame.length = sizeof(CalcResult);
			memcpy(frame.data, &calcResult, sizeof(CalcResult));
            //计算结果发送到客户端
			send_frame(socketClient, (char*)&frame, sizeof(frame));

			int sendedLength = 0;
            //将每一个三角面片的积分结果模值发送到客户端
			while (true)
			{
				strcpy(frame.command, "PreTriangleResult");
				if (sendedLength + 1024 < calcInfo.triangles_length*sizeof(float))
				{
					memcpy(frame.data, (char*)pre_triangle_result + sendedLength, 1024);
					frame.length = 1024;
					send_frame(socketClient, (char*)&frame, sizeof(Frame));
					sendedLength += 1024;
				}
				else
				{
					memcpy(frame.data, (char*)pre_triangle_result + sendedLength, calcInfo.triangles_length*sizeof(float) - sendedLength);
					frame.length = calcInfo.triangles_length*sizeof(float) - sendedLength;
					send_frame(socketClient, (char*)&frame, sizeof(Frame));
					break;
				}
			 }
			//将每一个三角面片的积分结果（实部虚部）发送到客户端
			sendedLength = 0;
			while (true)
			{
				strcpy(frame.command, "TriangleResultReIm");
				if (sendedLength + 1024 < 2*calcInfo.triangles_length*sizeof(float))
				{
					memcpy(frame.data, (char*)pre_triangle_result_reim + sendedLength, 1024);
					frame.length = 1024;
					send_frame(socketClient, (char*)&frame, sizeof(Frame));
					sendedLength += 1024;
				}
				else
				{
					memcpy(frame.data, (char*)pre_triangle_result_reim + sendedLength, 2*calcInfo.triangles_length*sizeof(float) - sendedLength);
					frame.length = 2*calcInfo.triangles_length*sizeof(float) - sendedLength;
					send_frame(socketClient, (char*)&frame, sizeof(Frame));
					break;
				}
			 }
			//将时域积分结果发送回客户端  姬梓遇20210913
			sendedLength = 0;
			while (true)
			{
				strcpy(frame.command, "re");
				if (sendedLength + 1024 < totalsize*sizeof(float))
				{
					memcpy(frame.data, (char*)result_1200 + sendedLength, 1024);
					frame.length = 1024;
					send_frame(socketClient, (char*)&frame, sizeof(Frame));
					sendedLength += 1024;
				}
				else
				{
					memcpy(frame.data, (char*)result_1200 + sendedLength, totalsize*sizeof(float) - sendedLength);
					frame.length = totalsize*sizeof(float) - sendedLength;
					send_frame(socketClient, (char*)&frame, sizeof(Frame));
					break;
				}
			}




			printf("%d\t%d\t%d\t%f\t%d\t%d\t%f\t%f\n", i, pre_device_width[i], height[i], calcInfo.config.pipe_size, 0, 0, result, time_all);
			//避免线程在其他地方cancel，以保证上位机能正常暂停计算 2022.3.24 jzy
			pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
			angle ++;          
			}
		
		}
		start_alpha = calcInfo.config.start_alpha;
		printf("第%d个阵元计算完毕\n",recv_index);
	}
    fclose(fileresult_TS);
	free(result_1200_re);
	free(result_1200_im);
	free(result_1200);
	//cudaFree(d_s_sum_re);
	//cudaFree(d_s_sum_re);
	//fclose(fileresult_1200);//20200919
	//cudaFree(plan.d_m_sum_re);
	//cudaFree(plan.d_m_sum_im);


	omp_set_num_threads(calcInfo.config.card_num);  // create as many CPU threads as there are CUDA devices
	#pragma omp parallel
	{
		int i = omp_get_thread_num();
		HANDLE_ERROR(cudaSetDevice(calcInfo.config.select_device_list[i]));

		FreeOnGPU(plan[i].d_rays1, plan[i].d_squares1, plan[i].d_rays2, plan[i].d_squares2, 
			plan[i].d_effrays, plan[i].d_center, plan[i].d_axis, plan[i].d_transMat, plan[i].d_reim, plan[i].d_sum_re, plan[i].d_sum_im, plan[i].d_sum_sre, plan[i].d_sum_sim,
			plan[i].d_DivRayTubeNum, plan[i].d_sum_gmem, plan[i].d_sum_Gmem, plan[i].d_squares_pred);

		free_data(plan[i].d_out_array, plan[i].d_tree, plan[i].d_points, plan[i].d_triangles);
	}

	//host Free

	if(plan!=NULL)
		free(plan);
	if(pre_device_height!=NULL)
		free(pre_device_height);
	if(pre_device_width!=NULL)
		free(pre_device_width);
	if(height!=NULL)
		free(height);
	if(e_fi_max!=NULL)
		free(e_fi_max);
	if(e_st_min!=NULL)
		free(e_st_min);
	if(h_out_array!=NULL)
		free(h_out_array);
	if(h_tree!=NULL)
		free(h_tree);
	
	if(h_reim!=NULL)
	{
		for(int i=0;i<calcInfo.config.card_num;i++)
			free(h_reim[i]);
		free(h_reim);
	}
	if(h_TSOfPerTriangle!=NULL)
	{
		for(int i=0;i<calcInfo.config.card_num;i++)
			free(h_TSOfPerTriangle[i]);
		free(h_TSOfPerTriangle);
	}
	if(pre_triangle_result!=NULL)
		free(pre_triangle_result);
	if(pre_triangle_result_reim!=NULL)
		free(pre_triangle_result_reim);
	if(SubAperturePlane!=NULL)
		free(SubAperturePlane);	
	if(AperturePlane!=NULL)
		free(AperturePlane);
	if(PreAngelTime!=NULL)
	{
		for(int i=0;i<2*calcInfo.config.card_num-1;i++)
			free(PreAngelTime[i]);
		free(PreAngelTime);
	}		
	if(New_receive_points!=NULL)
		free(New_receive_points);
	if(h_points!=NULL)
		free(h_points);
	Frame frame;
	strcpy(frame.command, "CalcOver");
	frame.length = 0;

	send_frame(socketClient, (char*)&frame, sizeof(frame));
	printf("calc over\n");
	pthread_cleanup_pop(0);


	return NULL;

}
