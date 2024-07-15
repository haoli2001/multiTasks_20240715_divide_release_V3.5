#ifndef COMMON_STRUCT
#define COMMON_STRUCT

#define DATA_Type float
#include <vector_types.h>
#define PI 3.14159265358979

/********************************
/**   ģ�ͽڵ���Ϣ
/*******************************/
struct Element
{
	DATA_Type point[3];//�±� 0 1 2 3�ֱ����õ�x,y,z��ֵ
	int PointsIndex;//�õ�ı��ֵ
};//16Bytes

/********************************
/**   ģ����Ԫ��Ϣ
/*******************************/
struct Triangle
{
	int Points[3];//�洢����ı��
	int TriangleIndex;//��������Ԫ�ı��
};//16Bytes

/********************************
/**   ��
/*******************************/
struct Axis//ֱ������
{
	float x;
	float y;
	float z;
};//12Bytes

struct Axis_slx//ֱ������
{
	float p[3];
};//12Bytes

/*******************************/
/**   ���߹��� ��������
/*******************************/
struct Square
{
	int4 CornerRayIndex; // ���߹����ĽǶ����߱��
	int CenterRayIndex; // ���߹������������߱��
	float3 CenterRay;
	bool right;
	bool IsDivRayTube;
};//49Bytes

/********************************
/**   ���߽ṹ
/*******************************/
struct Direction
{
	float3 p;           //������������Ԫ�Ľ���
	float3 dir;         //���ߵķ�������
	float distance; //�����뷢������
	int triangle_index;
	int times;
};

/********************************
/**   ��Χ��
/*******************************/
struct Box
{	//��Χ�е���������
	DATA_Type bmin[3];//������ԭ�����������ֵ x y z
	DATA_Type bmax[3];//������ԭ����Զ������ֵ x y z
};//24Bytes

/********************************
/**   ��Χ�� ��������
/*******************************/
struct Prim_Box
{
	DATA_Type bmin[3];//������ԭ�����������ֵ x y z
	DATA_Type bmax[3];//������ԭ����Զ������ֵ x y z
	DATA_Type bmid[3];//��Χ�е��е�
	int Box_Index; // ������Ԫ���
};//40Bytes


/********************************
/**   KDTree��״�ṹ�洢
/*******************************/
struct KD_Node                                
{
	int Split_Axis;//�ָ��� 0��ʾx,1��ʾy��2��ʾz
	int Depth;//��ǰ����
	int PrimCount;//������ԪƬ��
	int begin;
	int end;//���ýڵ�����������ʼ�±�ͽ����±�
	Element  SplitPos;//�ָ���λ��
	bool IsLeaf;//�Ƿ���Ҷ�ڵ�
	bool IsEmpty;
	struct Box box;//��Χ�е���Ϣ
	KD_Node *(rope[6]);//�������������Ϣ X������0��1�ң�Y������2��3�ң�Z������4��5��
	KD_Node *LeftChild;
	KD_Node *RightChild;

	int LeftIndex;
	int RightIndex;
	int RopeIndex[6];
	int index;
};

/********************************
/**  KD_Tree���Դ洢   
/*******************************/
struct KD_Node_V                                 
{
	int Split_Axis;//�ָ��� 0��ʾx,1��ʾy��2��ʾz
	int PrimCount;//������ԪƬ��
	int begin;
	int end;//���ýڵ�����������ʼ�±�ͽ����±�
	Element  SplitPos;//�ָ���λ��
	bool IsLeaf;//�Ƿ���Ҷ�ڵ�
	bool IsEmpty;
	struct Box box;
	int LeftIndex;
	int RightIndex;
	int RopeIndex[6];
	int index;
};//94Bytes


/********************************
/**  ����ͨ��֡�ṹ
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
  float sampling_width;//��������
  float tao;
  float time_start_frequency;
  float time_end_frequency;
  float relative_velocity;
  float integral_gain;
  float velocity1;
  float velocity2;
  float reflect_coeff;
  bool reflect_coeff_Auto_flag; //����ϵ������ģʽ�� tureΪ�Զ����㣬 falseΪ�û�����ֵ
  float H_axis_angle;			//��ͨ����λˮƽ���
	float V_axis_angle;			//��ͨ����λ��ֱ���
  int recvPointsNum;          //��ͨ�����յ���Ŀ
  int recvPointsStartIdx;      //��ͨ�����յ������ʼ���
  float continue_alpha;		//��λ����ͣ�����󣬼������е���ʼalphaֵ(-1��ʾ��������)
  float water_line;    //ˮ�����
	float sway_theta;    //����ζ��������Ƕ�
	float sway_phi;
	float sway_psi;
};

struct CalcResult
{
	float angle;                       //�˽���ĽǶ�
	float freq;
	
	float TS;          //����TSǿ��ֵ
	
	float calc_time;
	int raysnum;                     //10��28������
	int squarenum;					 //10��28������
	int height;
	int width;
  int recvIdx;             //��¼��ǰ����Ľ��յ���� 22.9.5
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
