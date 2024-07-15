#ifndef COMMON_STRUCT
#define COMMON_STRUCT

#define DATA_Type float
#include <vector_types.h>
#define PI 3.14159265358979

/********************************
/**   模型节点信息
/*******************************/
struct Element
{
	DATA_Type point[3];//下标 0 1 2 3分别代表该点x,y,z的值
	int PointsIndex;//该点的编号值
};//16Bytes

/********************************
/**   模型面元信息
/*******************************/
struct Triangle
{
	int Points[3];//存储顶点的编号
	int TriangleIndex;//该三角面元的编号
};//16Bytes

/********************************
/**   点
/*******************************/
struct Axis//直角坐标
{
	float x;
	float y;
	float z;
};//12Bytes

struct Axis_slx//直角坐标
{
	float p[3];
};//12Bytes

/*******************************/
/**   声线管束 （索引）
/*******************************/
struct Square
{
	int4 CornerRayIndex; // 声线管束的角顶射线编号
	int CenterRayIndex; // 声线管束的中心射线编号
	float3 CenterRay;
	bool right;
	bool IsDivRayTube;
};//49Bytes

/********************************
/**   射线结构
/*******************************/
struct Direction
{
	float3 p;           //射线与三角面元的交点
	float3 dir;         //射线的方向向量
	float distance; //交点与发射点距离
	int triangle_index;
	int times;
};

/********************************
/**   包围盒
/*******************************/
struct Box
{	//包围盒的两个顶点
	DATA_Type bmin[3];//离坐标原点最近的三个值 x y z
	DATA_Type bmax[3];//离坐标原点最远的三个值 x y z
};//24Bytes

/********************************
/**   包围盒 （进化）
/*******************************/
struct Prim_Box
{
	DATA_Type bmin[3];//离坐标原点最近的三个值 x y z
	DATA_Type bmax[3];//离坐标原点最远的三个值 x y z
	DATA_Type bmid[3];//包围盒的中点
	int Box_Index; // 三角面元序号
};//40Bytes


/********************************
/**   KDTree数状结构存储
/*******************************/
struct KD_Node                                
{
	int Split_Axis;//分割轴 0表示x,1表示y，2表示z
	int Depth;//当前树高
	int PrimCount;//三角面元片数
	int begin;
	int end;//给该节点分配的数组起始下标和结束下标
	Element  SplitPos;//分割点的位置
	bool IsLeaf;//是否是叶节点
	bool IsEmpty;
	struct Box box;//包围盒的信息
	KD_Node *(rope[6]);//六个面的线索信息 X方向上0左1右，Y方向上2左3右，Z方向上4左5右
	KD_Node *LeftChild;
	KD_Node *RightChild;

	int LeftIndex;
	int RightIndex;
	int RopeIndex[6];
	int index;
};

/********************************
/**  KD_Tree线性存储   
/*******************************/
struct KD_Node_V                                 
{
	int Split_Axis;//分割轴 0表示x,1表示y，2表示z
	int PrimCount;//三角面元片数
	int begin;
	int end;//给该节点分配的数组起始下标和结束下标
	Element  SplitPos;//分割点的位置
	bool IsLeaf;//是否是叶节点
	bool IsEmpty;
	struct Box box;
	int LeftIndex;
	int RightIndex;
	int RopeIndex[6];
	int index;
};//94Bytes


/********************************
/**  网络通信帧结构
/*******************************/
struct Frame
{
	char command[20];
	int length;
	char data[1024];
};

struct ConfigStruct
{
	float start_beta;
	float start_alpha;
	float end_alpha;
	float wave_length;
	float pipe_size;
	float start_frequency;
	float end_frequency;
	float far_distance;
	int card_num;
	int select_device_list[20];
  float sampling_frequency;
  float sampling_width;//采样长度
  float tao;
  float time_start_frequency;
  float time_end_frequency;
  float relative_velocity;
  float integral_gain;
  float velocity1;
  float velocity2;
  float reflect_coeff;
  bool reflect_coeff_Auto_flag; //反射系数配置模式， ture为自动计算， false为用户输入值
  float H_axis_angle;			//多通道阵位水平倾角
	float V_axis_angle;			//多通道阵位垂直倾角
  int recvPointsNum;          //多通道接收点数目
  int recvPointsStartIdx;      //多通道接收点计算起始序号
  float continue_alpha;		//上位机暂停重启后，继续运行的起始alpha值(-1表示正常运行)
  float water_line;    //水线深度
	float sway_theta;    //船体晃动的三个角度
	float sway_phi;
	float sway_psi;
};

struct CalcResult
{
	float angle;                       //此结果的角度
	float freq;
	
	float TS;          //积分TS强度值
	
	float calc_time;
	int raysnum;                     //10月28，孙力
	int squarenum;					 //10月28，孙力
	int height;
	int width;
  int recvIdx;             //记录当前计算的接收点序号 22.9.5
};

struct GPUWatchStruct
{
	int device_id;
	int gpu;
	int memory;
	int temp;
	int shutdown_temp;
	int slowdown_temp;
	long long int total;
	long long int used;
	long long int free;
};

struct DeviceInfo
{
	int deviceID;
	char deviceName[30];
	int deviceCount;
	int coresPreMutiprocess;
	int mutiprocessCount;

};
#endif
