#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include "kd_struct.h"

int Mid_Choose_Split_Axis(KD_Node *treenode, Prim_Box *array, int arr_length, Element* mid);
bool Cut_Blank(KD_Node *treenode, Prim_Box **array, Prim_Box **arr_left, Prim_Box **arr_right, int arr_length, int *nl, int *nr);
void Count_Prim(Prim_Box *array, int arr_length, DATA_Type Split_Pos, int *nl, int *nr, int *np, int *Axis);
void Count_Prim_Mid(Prim_Box *array, int arr_length, DATA_Type Split_Pos, int *nl, int *nr, int Axis);
int SAH_Choose_Split_Axis(KD_Node *treenode, Prim_Box *array, int arr_length,int *box_flag, DATA_Type *Position,int *np);
void Updata_Node(KD_Node *treenode, Prim_Box *array);
void StopBuild(KD_Node * treenode, Prim_Box *array, Prim_Box *out_array, int arr_length, int *out_arr_length);
//void Preprocessing_Triangles(FILE *fp, Triangle *Triangles, Element *Points, struct Prim_Box *Prim, int Triangle_Num);
void Max_Min(Prim_Box *array, int arr_length, int axis, DATA_Type *max, DATA_Type *min);
void Resize_Array(Prim_Box *array, Prim_Box *arr_left, Prim_Box *arr_right, int arr_length, int nl, int nr, int flag, DATA_Type Split_Pos, int Axis);
void Resize_Array(Prim_Box *array, Prim_Box *arr_left, Prim_Box *arr_right, int arr_length, int nl, int nr, DATA_Type Split_Pos, int Axis);
void Resize_Array_Mid(Prim_Box *array, Prim_Box *arr_left, Prim_Box *arr_right, int arr_length, int nl, int nr, DATA_Type Split_Pos, int Axis);


/********************************************************************************************************
函数功能：对三角面元进行预处理，每个面元都使用包围盒包起来，便于SAH计算
参数含义：1、要打开的文件指针
2、建立的三角面元的数组
3、保存顶点信息
4、保存建立的面元包围盒
********************************************************************************************************/
void Preprocessing_Triangles(FILE *fp, Triangle *Triangles, Element *Points, int *Node_Num,int *Triangle_Num)
{
	int i, j;
	int temp;//保存后两位无用数据

/*-----读取dat文件------*/
	int node_num;
	int face_num;
	int index = 0,num = 0;
	int tempnum1,tempnum2,tempnum3;

	fscanf(fp,"%d",&num);
	printf("node_num=%d\n",num);						// 读节点数;
	float *pointx = new float[num];					// 节点的x坐标
	float *pointy = new float[num];					// 节点的y坐标
	float *pointz = new float[num];					// 节点的z坐标
	while (num)
	{

	        fscanf(fp,"%d",&index);

		///////球体和柱体适用/////
		fscanf(fp,"%f",&(pointx[index-1]));			// 读x坐标;     
		fscanf(fp,"%f",&(pointy[index-1]));			// 读y坐标;
		fscanf(fp,"%f",&(pointz[index-1]));			// 读z坐标;
                //////////////////////////
		
		///////仅(benchmark)适用/////
		//fscanf(fp,"%d",&index);
		//fscanf(fp,"%f",&(pointz[index-1]));			// 读x坐标;     
		//fscanf(fp,"%f",&(pointy[index-1]));			// 读y坐标;
		//fscanf(fp,"%f",&(pointx[index-1]));			// 读z坐标;
		//pointx[index-1] = pointx[index-1]  + 24.0*1.5;
		//pointx[index-1] = pointx[index-1] ;
		//pointy[index-1] = pointy[index-1] ;
		//pointz[index-1] = pointz[index-1] ;
		//////////////////////////////

		num--;
	}
	node_num = index;                                               //节点个数

	fscanf(fp,"%d",&num);	
	printf("ele_num=%d\n",num);					// 读面元数;
	int *faceOne = new int[num];					// 节点的x坐标
	int *faceTwo = new int[num];					// 节点的y坐标
	int *faceThree = new int[num];					// 节点的z坐标
	while (num)
	{
		fscanf(fp,"%d",&index);
		fscanf(fp,"%d",&(faceOne[index-1]));			// 读x坐标;     
		fscanf(fp,"%d",&(faceTwo[index-1]));			// 读y坐标;
		fscanf(fp,"%d",&(faceThree[index-1]));			// 读z坐标;
		num--;
	}
	
	face_num = index;                                                //面元个数




	for (i = 0; i < node_num; i++)
	{
		
			/*(Points + i)->point[0] = pointx[i];
			(Points + i)->point[1] = pointy[i];
			(Points + i)->point[2] = pointz[i];*/

	        Points[i].point[0] = pointx[i];
			Points[i].point[1] = pointy[i];
			Points[i].point[2] = pointz[i];

		//fscanf(fp,"%d",&(Points + i)->PointsIndex);//给顶点进行编号
		//(Points + i)->PointsIndex = i;
                Points[i].PointsIndex = i;
		//printf("\nRead  Point%d", (Points + i)->PointsIndex);
	}

	for (i = 0; i<face_num; i++)
	{

		/*(Triangles + i)->Points[0] = faceOne[i] ;
		(Triangles + i)->Points[1] = faceTwo[i] ;
		(Triangles + i)->Points[2] = faceThree[i] ;

		(Triangles + i)->TriangleIndex = i;*/

        Triangles[i].Points[0] = faceOne[i] - 1;
		Triangles[i].Points[1] = faceTwo[i] - 1;
		Triangles[i].Points[2] = faceThree[i] - 1;

		Triangles[i].TriangleIndex = i;


	}
	*Node_Num = node_num;
	*Triangle_Num = face_num;
	//Sort_Box(Points, Triangles, Prim, Triangle_Num);
}

/*void Preprocessing_Triangles(FILE *fp, Triangle *Triangles, Element *Points, Prim_Box *Prim, int Triangle_Num)
{
	int i, j;
	int temp;//保存后两位无用数据


	for (i = 0; i < Points_Num; i++)
	{
		for (j = 0; j < 3; j++)
			fscanf(fp, "%f", &(Points + i)->point[j]);
		//fscanf(fp,"%d",&(Points + i)->PointsIndex);//给顶点进行编号
		(Points + i)->PointsIndex = i;
		//printf("\nRead  Point%d", (Points + i)->PointsIndex);
	}

	for (i = 0; i<Triangle_Num; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(fp,"%x",&(Triangles + i)->Points[j]);//三角面元内存储的顶点编号
			(Triangles + i)->Points[j] = (Triangles + i)->Points[j] -1 ;//由于顶点索引值由1开始，因此要减一
		}
		fscanf(fp, "%d", &temp);//三角面元内存储的顶点编号
		fscanf(fp, "%d", &temp);//三角面元内存储的顶点编号
		(Triangles + i)->TriangleIndex = i;
		//fscanf(fp,"%d",&(Triangles + i)->TriangleIndex);//三角形编号
		//printf("\nRead  Triangle%d", (Triangles + i)->TriangleIndex);
	}

	Sort_Box(Points, Triangles, Prim, Triangle_Num);
}*/
/*******************************************************
函数功能：算出所有三角面元的包围盒的两个顶点的坐标值
参数含义:1、包含所有三角面元的数组
2、包含所有顶点的数组
3、由所有面元包围盒构成的数组
******************************************************/
void Sort_Box(Element points[], Triangle triangle[], Prim_Box *prim, int Triangle_Num)
{

	int i, j;
	DATA_Type Max, Min;
	for (i = 0; i < Triangle_Num; i++)
	{
		for (j = 0; j < 3; j++)
		{
			Max = (points + (triangle + i)->Points[0])->point[j] >(points + (triangle + i)->Points[1])->point[j] ? (points + (triangle + i)->Points[0])->point[j] : (points + (triangle + i)->Points[1])->point[j];
			Max = Max >(points + (triangle + i)->Points[2])->point[j] ? Max : (points + (triangle + i)->Points[2])->point[j];
			Min = (points + (triangle + i)->Points[0])->point[j] < (points + (triangle + i)->Points[1])->point[j] ? (points + (triangle + i)->Points[0])->point[j] : (points + (triangle + i)->Points[1])->point[j];
			Min = Min < (points + (triangle + i)->Points[2])->point[j] ? Min : (points + (triangle + i)->Points[2])->point[j];
			//if (Max != Min)
			//{
				(prim + i)->bmax[j] = Max;
				(prim + i)->bmin[j] = Min;
			//}
			//else//当三角形垂直于坐标轴时，给其一个宽度
			//{
			//	(prim + i)->bmax[j] = Max+Lambda/2;
			//	(prim + i)->bmin[j] = Min - Lambda/2;
			//}
		}
		(prim + i)->Box_Index = i;
		//每个包围盒中三角形的中心
		(prim + i)->bmid[0] = ((points + (triangle + i)->Points[0])->point[0] + (points + (triangle + i)->Points[1])->point[0] + (points + (triangle + i)->Points[2])->point[0]) / 3;
		(prim + i)->bmid[1] = ((points + (triangle + i)->Points[0])->point[1] + (points + (triangle + i)->Points[1])->point[1] + (points + (triangle + i)->Points[2])->point[1]) / 3;
		(prim + i)->bmid[2] = ((points + (triangle + i)->Points[0])->point[2] + (points + (triangle + i)->Points[1])->point[2] + (points + (triangle + i)->Points[2])->point[2]) / 3;
	}

}
/****************************************************************************************
函数功能：对节点进行初始化
参数意义：1、初始节点
2、所有面元包围盒构成的数组
返回值意义：选择的分割轴 0.1.2分别代表x.y.z轴
****************************************************************************************/
void KD_Node_init(struct KD_Node *kd_node, Prim_Box *array,int arr_length)
{//节点初始化
	//kd_node = (struct KD_Node*)(malloc( sizeof( struct KD_Node ) ));
	kd_node->Split_Axis = 0;
	kd_node->Depth = 1;
	kd_node->PrimCount = arr_length;
	kd_node->begin = 0;
	kd_node->end = arr_length - 1;
	//kd_node->NodeIndex = 0;
	for (int i = 0; i < 3; i++)
		kd_node->SplitPos.point[i] = 0;
	kd_node->SplitPos.PointsIndex = 0;
	kd_node->IsLeaf = false;
	kd_node->IsEmpty = false;
	for (int i = 0; i<6; i++)
		kd_node->rope[i] = NULL;
	//给初始节点的包围盒确定范围，即找到包围盒中离原点最近和最远的点
	Max_Min(array, arr_length, 0, &kd_node->box.bmax[0], &kd_node->box.bmin[0]);
	Max_Min(array, arr_length, 1, &kd_node->box.bmax[1], &kd_node->box.bmin[1]);
	Max_Min(array, arr_length, 2, &kd_node->box.bmax[2], &kd_node->box.bmin[2]);
	return;
}
/******************************** 
 *函数名：swap 
 *作用：交换两个结构体的值 
 *参数：交换的两个结构体 
 *返回值：无 
 ********************************/  
void swap(Prim_Box *a, Prim_Box *b)    
{  
    Prim_Box temp;  
    temp = *a;  
    *a = *b;  
    *b = temp;  
    return ;  
}  
  
/************************************ 
 *函数名：quicksort 
 *作用：快速排序算法，从小到大排序
 *参数： 需排序的结构体数组 ，数组中的起始位置，数组中结束位置，需要排序的维度
 *返回值：无 
 ************************************/  
void QuickSort(Prim_Box* array,  int begin, int end,int dim) 
{  
    int i, j;  
    if(begin < end)  
    {  
        i = begin + 1;  // 将array[begin]作为基准数，因此从array[begin+1]开始与基准数比较！  
        j = end;        // array[end]是数组的最后一位  
            
        while(i < j)  
        {  
            if(  (array+i)->bmid[dim] > (array+begin) ->bmid[dim])  // 如果比较的数组元素大于基准数，则交换位置。  
            {  
                swap(&array[i],&array[j]);  // 交换两个结构体
                j--;  
            }  
            else  
            {  
                i++;  // 将数组向后移一位，继续与基准数比较。  
            }  
        }  
		/* 跳出while循环后，i = j。 
         * 此时数组被分割成两个部分  -->  array[begin+1] ~ array[i-1] < array[begin] 
         *                           -->  array[i+1] ~ array[end] > array[begin] 
         * 这个时候将数组array分成两个部分，再将array[i]与array[begin]进行比较，决定array[i]的位置。 
         * 最后将array[i]与array[begin]交换，进行两个分割部分的排序！以此类推，直到最后i = j不满足条件就退出！ 
         */  
        if( (array+i)->bmid[dim] >= (array+begin)->bmid[dim])  // 这里必须要取等“>=”，否则数组元素由相同的值时，会出现错误！  
        {  
            i--;  
        }  
        swap(&array[begin], &array[i]);  // 交换array[i]与array[begin]  
        QuickSort(array, begin, i,dim);  
        QuickSort(array, j, end,dim);  
    }  
}  

/****************************************************************************************
函数功能：使用中分法选择分割面
参数意义：1、由包围盒构成的数组
2、3、begin和end表示本次排序在数组中的起始点
4、mid表示分割点的值
返回值意义：选择的分割轴 0.1.2分别代表x.y.z轴
****************************************************************************************/
int Mid_Choose_Split_Axis(KD_Node *treenode, Prim_Box *array, int arr_length, Element* mid)
{
	//三个方向上的长度
	DATA_Type x_dim = 0;
	DATA_Type y_dim = 0;
	DATA_Type z_dim = 0;
	//三个方向上的中间值，即分割点
	DATA_Type x_mid = 0;
	DATA_Type y_mid = 0;
	DATA_Type z_mid = 0;
	int nl,nr;
	//三个方向上排序后的数组
	//int x_y_z = 0;//0表示x方向，1表示y，2表示z
	/***********************************************************
	求三个维度上的长度，再找出最长的那个方向，并记录中点值
	***********************************************************/

	x_dim = fabs(treenode->box.bmax[0] - treenode->box.bmin[0]);
	x_mid = x_dim / 2 + treenode->box.bmin[0];

	y_dim = fabs(treenode->box.bmax[1] - treenode->box.bmin[1]);
	y_mid = y_dim / 2 + treenode->box.bmin[1];

	z_dim = fabs(treenode->box.bmax[2] - treenode->box.bmin[2]);
	z_mid = z_dim / 2 + treenode->box.bmin[2];

	if (x_dim >= y_dim && x_dim >= z_dim)
	{
		mid->point[0] = x_mid;
		//Count_Prim_Mid(array, begin, end, x_mid, &nl, &nr, 0);
		//treenode->SplitPos.PointsIndex = nl;
		//Resize_Array_Mid(array, begin, end, nl, nr, x_mid, 0);
		return 0;
	}
	else if (y_dim >= x_dim && y_dim >= z_dim)
	{
		mid->point[1] = y_mid;
		//Count_Prim_Mid(array, begin, end, y_mid, &nl, &nr, 1);
		//treenode->SplitPos.PointsIndex = nl;
		//Resize_Array_Mid(array, begin, end, nl, nr, y_mid, 1);
		return 1;
	}
	else if (z_dim >= y_dim && z_dim >= x_dim)
	{
		mid->point[2] = z_mid;
		//Count_Prim_Mid(array, begin, end, z_mid, &nl, &nr, 2);
		//treenode->SplitPos.PointsIndex = nl;
		//Resize_Array_Mid(array, begin, end, nl, nr, z_mid, 2);
		return 2;
	}
	return -1;
}
/*********************************************************
函数功能：剔除包围盒中的空白部分

*********************************************************/
bool Cut_Blank(KD_Node *treenode, Prim_Box **array, Prim_Box **arr_left, Prim_Box **arr_right, int arr_length, int *nl, int *nr)
{
	DATA_Type max, min, blank, length,empty_rate;
	empty_rate = 0.20;
	for (int i = 0; i < 3; i++)//三个维度上找空白部分，当超过阈值时进行分割
	{
		Max_Min(*array, arr_length, i, &max, &min);
		length = max - min;
		blank = fabs(max - treenode->box.bmax[i]);//x轴上右孩子为空节点
		if (blank / length > empty_rate)
		{
			treenode->RightChild->IsEmpty = true;
			treenode->RightChild->PrimCount = 0;

			treenode->LeftChild->PrimCount = arr_length;
			//arr_left = (Prim_Box*)malloc(arr_length*sizeof(Prim_Box));
			*arr_left = *array;//右孩子为空，将当前数组直接转给左孩子
			*nl = arr_length;
			*nr = 0;
			treenode->SplitPos.point[i] = max;
			treenode->Split_Axis = i;
			Updata_Node(treenode, *array);
			return true;
		}
		blank = fabs(min - treenode->box.bmin[i]);
		if (blank / length > empty_rate)
		{
			treenode->LeftChild->IsEmpty = true;
			treenode->LeftChild->PrimCount = 0;
			//treenode->RightChild->begin = treenode->begin;
			//treenode->RightChild->end = treenode->end;
			treenode->RightChild->PrimCount = arr_length;
			//arr_right = (Prim_Box*)malloc(arr_length*sizeof(Prim_Box));
			*arr_right = *array;
			*nl = 0;
			*nr = arr_length;
			treenode->SplitPos.point[i] = min;
			treenode->Split_Axis = i;
			Updata_Node(treenode, *array);
			return true;
		}

	}
	return false;
}
/*******************************************************
函数功能：找出数组中的最大最小值
参数意义：
1、数组名
2、数组中起点
3、待查找分割轴
4、最大值
5、最小值
********************************************************/
void Max_Min(Prim_Box *array, int arr_length, int axis, DATA_Type *max, DATA_Type *min)
{

	*min = array ->bmin[axis];
	*max =array ->bmax[axis];
	for (int i = 0; i < arr_length; i++)
	{
		if ((array + i)->bmin[axis] < *min)
			*min = (array + i)->bmin[axis];
		if ((array + i)->bmax[axis] > *max)
			*max = (array + i)->bmax[axis];
	}
}
/****************************************************************************************
函数功能：输出分割平面左右以及被分割平面穿过的子包围盒的数量
参数意义：1、由包围盒构成的数组
2、3、begin和end表示本次排序在数组中的起始点
4、Poisition内包含了分割点位置以及分割平面所在包围盒的编号
5、该节点中包围盒最大值小于分割点的面元数量
6、该节点中包围盒最大值大于分割点的面元数量
7、该节点中包围盒被分割平面分割的面元数量
8、分割轴
返回值意义：无返回值
****************************************************************************************/
void Count_Prim(Prim_Box *array, int arr_length, DATA_Type Split_Pos, int *nl, int *nr, int *np, int *Axis)
{
	*nr = 0;
	*nl = 0;
	*np = 0;

	for (Prim_Box *itr = array; itr != (array + arr_length); itr++)
	{
		if (itr->bmax[*Axis] == itr->bmin[*Axis] && (itr->bmax[*Axis]) == Split_Pos)
		{
			//(*nl)++;
			//(*nr)++;
			continue;
		}
		else
		{
			(*nl) += ((itr->bmax[*Axis]) <= Split_Pos ? 1 : 0);//待修改
			(*nr) += ((itr->bmin[*Axis]) >= Split_Pos ? 1 : 0);
		}

	}
	*np = arr_length - *nr - *nl ;
	return;
}
/****************************************************************************************
函数功能：输出分割平面左右以及被分割平面穿过的子包围盒的数量
参数意义：1、由包围盒构成的数组
2、3、begin和end表示本次排序在数组中的起始点
4、Poisition内包含了分割点位置以及分割平面所在包围盒的编号
5、该节点中包围盒最大值小于分割点的面元数量
6、该节点中包围盒最大值大于分割点的面元数量
7、该节点中包围盒被分割平面分割的面元数量
8、分割轴
返回值意义：无返回值
****************************************************************************************/
void Count_Prim_Mid(Prim_Box *array, int arr_length, DATA_Type Split_Pos, int *nl, int *nr, int Axis)
{
	*nr = 0;
	*nl = 0;
	for (Prim_Box *itr = array ; itr != (array + arr_length); itr++)
	{
		(*nl) += ((itr->bmin[Axis]) <= Split_Pos ? 1 : 0);//待修改
		(*nr) += ((itr->bmax[Axis]) >= Split_Pos ? 1 : 0);
	}
	return;
}
/****************************************************************************************
函数功能：将数组按照剖分平面左右切分，并覆盖原来的位置
参数含义：
1、数组名
2、3、待分离的面元在数组中的起始位置与终止位置
4、5、分别属于左右节点面元的数量
6、重叠面元归属左边还是右边的标志位 0是左边 1是右边
7、分割位置
8、分割轴
****************************************************************************************/
void Resize_Array(Prim_Box *array, Prim_Box *arr_left, Prim_Box *arr_right, int arr_length, int nl, int nr, int flag, DATA_Type Split_Pos, int Axis)
{
	//struct Prim_Box *left = (Prim_Box*)malloc(nl*(sizeof(Prim_Box)));
	//struct Prim_Box *right = (Prim_Box*)malloc(nr*(sizeof(Prim_Box)));
	//memset(left, 0, nl*sizeof(struct Prim_Box));//初始化为0
	//memset(right, 0, nr*sizeof(struct Prim_Box));//初始化为0

	if (!flag)//重叠面元在左边
		for (int i = 0; i< arr_length; i++)
		{
			if (((array + i)->bmin[Axis]) >= Split_Pos)
			{
				*arr_right = *(array + i);
				arr_right++;
			}
			else
			{
				*arr_left = *(array + i);
				arr_left++;
			}
		}
	else
		for (int i = 0; i < arr_length; i++)
		{
			if ((array + i)->bmax[Axis] <= Split_Pos)
			{
				*arr_left = *(array + i);
				arr_left++;
			}
			else
			{
				*arr_right = *(array + i);
				arr_right++;
			}
		}
	//array = NULL;
	//free(array);
}
void Resize_Array(Prim_Box *array, Prim_Box *arr_left, Prim_Box *arr_right, int arr_length, int nl, int nr, DATA_Type Split_Pos, int Axis)
{
	//struct Prim_Box *left = (Prim_Box*)malloc(nl*(sizeof(Prim_Box)));
	//struct Prim_Box *right = (Prim_Box*)malloc(nr*(sizeof(Prim_Box)));
	//memset(left, 0, nl*sizeof(struct Prim_Box));//初始化为0
	//memset(right, 0, nr*sizeof(struct Prim_Box));//初始化为0
	int l_index, r_index;
	l_index = 0;
	r_index = 0;
	for (int i = 0; i< arr_length; i++)
	{
		if ((array[i].bmax[Axis]) == array[i].bmin[Axis] && (array[i].bmax[Axis]) == Split_Pos)
		{
			arr_right[r_index] = array[i];
			r_index++;
			arr_left[l_index] = array[i];
			l_index++;
			continue;
		}
		else if ((array[i].bmax[Axis]) <= Split_Pos)
		{
			arr_left[l_index] = array[i];
			l_index++;
		}
		else if (array[i].bmin[Axis] >= Split_Pos)
		{
			arr_right[r_index] = array[i];
			r_index++;
		}
		else
		{
			arr_right[r_index] = array[i];
			r_index++;
			arr_left[l_index] = array[i];
			l_index++;
		}
	}
	//printf(" ");
}
void Resize_Array_Mid(Prim_Box *array, Prim_Box *arr_left, Prim_Box *arr_right, int arr_length, int nl, int nr,DATA_Type Split_Pos, int Axis)
{
	//struct Prim_Box *left = (Prim_Box*)malloc(nl*(sizeof(Prim_Box)));
	//struct Prim_Box *right = (Prim_Box*)malloc(nr*(sizeof(Prim_Box)));
	//memset(left, 0, nl*sizeof(struct Prim_Box));//初始化为0
	//memset(right, 0, nr*sizeof(struct Prim_Box));//初始化为0
	int l_index, r_index;
	l_index = 0;
	r_index = 0;
	for (int i = 0; i< arr_length; i++)
		{
			if ((array[i].bmax[Axis]) >= Split_Pos)
			{
				arr_right[r_index] = array [i];
				r_index++;
			}
			if (array[i].bmin[Axis]<= Split_Pos)
			{
				arr_left[l_index] = array[i];
				l_index++;
			}
		}
}
/****************************************************************************************
函数功能：使用SAH算法选择分割面
参数意义：1、由包围盒构成的数组
2、3、begin和end表示本次排序在数组中的起始点
4、Poisition内包含了分割点位置以及分割平面所在包围盒的编号
返回值意义：选择的分割轴 0.1.2分别代表x.y.z轴
****************************************************************************************/
int SAH_Choose_Split_Axis(KD_Node *treenode, Prim_Box *array, int arr_length, int *box_flag,struct Element *Position,int *np)
{
	int flag = -1;//选取的坐标轴
	//int box_flag = 0;
	int temp_flag = 0;
	DATA_Type temp1, temp2;//选择节点左右侧的标志
	int c_hit, c_walk,NL,NR,NP;
	c_hit = 19;
	c_walk = 9;//??????????????????????????????????????????
	//int alpha_k = 5;//alpha_k=c_hit/c_walk;
	DATA_Type cost = 0;
	DATA_Type min_cost, min_cost_L, min_cost_R;//选择的剖分面是选中包围盒的左边还是右边
	DATA_Type Surface,Surface_L,Surface_R;
	DATA_Type H, L, W,
					   H_L,L_L,W_L,
					   H_R,L_R,W_R;
	DATA_Type max, min;
		//当前节点包围盒的长宽高
		//printf("\nStart Split");
	//Max_Min(array, begin, end, 0, &max, &min);
	L = DataAbs(treenode->box.bmax[0] - treenode->box.bmin[0]);
	//treenode->box.bmin[0] = min;
	//treenode->box.bmax[0] = max;
	//Max_Min(array, begin, end, 1, &max, &min);
	W = DataAbs(treenode->box.bmax[1] - treenode->box.bmin[1]);
	//treenode->box.bmin[1] = min;
	//treenode->box.bmax[1] = max;
	//Max_Min(array, begin, end, 2, &max, &min);
	H = DataAbs(treenode->box.bmax[2] - treenode->box.bmin[2]);
	//treenode->box.bmin[2] = min;
	//treenode->box.bmax[2] = max;
	Surface = (H*L+H*W+W*L);
	min_cost = DATA_Type(c_hit * 10 * (arr_length)+c_walk);//初始化最小的代价为最大值
	for (int i = 0; i < 3; i++)
	{
		min_cost_L = min_cost;
		min_cost_R = min_cost;
		//QuickSort(array, begin, end, i);
		for (int j = 0; j <arr_length; j++)
		{
			if ((array + j)->bmin[i] < treenode->box.bmin[i] || (array + j)->bmax[i] > treenode->box.bmax[i])
				continue;
			switch (i)//以包围盒的左平面作为剖分面，被选做剖分节点的包围盒被划到右节点中
			{
				case 0:
				{
					L_L = DataAbs((array + j)->bmin[0] - treenode->box.bmin[0]);
					L_R = DataAbs(treenode->box.bmax[0] - (array + j)->bmin[0]);
					Surface_L =  (H*W + H*L_L + L_L*W);
					Surface_R =  (H*W + H*L_R + L_R*W);
					break;
				}
				case 1:
				{
					W_L = DataAbs((array + j)->bmin[1] - treenode->box.bmin[1]);
					W_R = DataAbs(treenode->box.bmax[1] - (array + j)->bmin[1]);
					Surface_L =  (H*W_L + H*L + L*W_L);
					Surface_R =  (H*W_R + H*L + L* W_R); 
					break;
				}
				case 2:
				{
					H_L = DataAbs((array + j)->bmin[2] - treenode->box.bmin[2]);
					H_R = DataAbs(treenode->box.bmax[2] - (array + j)->bmin[2]);
					Surface_L =  (H_L*W + H_L*L + L*W);
					Surface_R =  (H_R*W + H_R*L + L*W);
					break;
				}
			}
			//printf("\nStart count prim:%d",j);
			Count_Prim(array, arr_length, (array + j)->bmin[i], &NL, &NR, &NP, &i);
			temp1 = (c_walk + (Surface_L*(NL + NP)*c_hit + Surface_R*NR*c_hit) / Surface);
			temp2 = (c_walk + (Surface_L*NL*c_hit + Surface_R*(NR + NP)*c_hit) / Surface);
			if (temp1 <= temp2 )//重叠的包围盒是划到左边还是右边
			{
				min_cost_L = temp1;
				//NL += NP;
				temp_flag = 0;//重叠包围盒划分标志
			}
			else
			{
				min_cost_L = temp2;
				temp_flag = 1;
			}
			if (min_cost_L < min_cost)
			{
				Position->point[i] = (array + j)->bmin[i];//记录分割点位
				Position->PointsIndex = NL ;//若重叠节点归类到左边，则包围盒大小与原来不一致，边界变大
				*np = NP;
				min_cost = min_cost_L;
				flag = i;//当在某维度上出现更小的代价时，记录其维度
				*box_flag = temp_flag;
			}
			//min_cost_L = min_cost_L <= cost ? min_cost_L : cost;
			switch (i)//以包围盒的右平面作为剖分面，被选做剖分节点的包围盒被划到左节点中
			{
				case 0:
				{
					L_L = DataAbs((array + j)->bmax[0] - treenode->box.bmin[0]);
					L_R = DataAbs(treenode->box.bmax[0] - (array + j)->bmax[0]);
					Surface_L =  (H*W + H*L_L + L_L*W);
					Surface_R =  (H*W + H*L_R + L_R*W);
					break;
				}
				case 1:
				{
					W_L = DataAbs((array + j)->bmax[1] - treenode->box.bmin[1]);
					W_R = DataAbs(treenode->box.bmax[1] - (array + j)->bmax[1]);
					Surface_L =  (H*W_L + H*L + L*W_L);
					Surface_R =  (H*W_R + H*L + L* W_R);
					break;
				} 
				case 2:
				{
					H_L = DataAbs((array + j)->bmax[2] - treenode->box.bmin[2]);
					H_R = DataAbs(treenode->box.bmax[2] - (array + j)->bmax[2]);
					Surface_L =  (H_L*W + H_L*L + L*W);
					Surface_R =  (H_R*W + H_R*L + L*W);
					break;
				}
			}
			Count_Prim(array, arr_length, (array + j)->bmax[i], &NL, &NR, &NP, &i);
			temp1 = (c_walk + (Surface_L*(NL + NP)*c_hit + Surface_R*NR*c_hit) / Surface);
			temp2 = (c_walk + (Surface_L*NL*c_hit + Surface_R*(NR + NP)*c_hit) / Surface);
			if (temp1 <= temp2)//重叠的包围盒是划到左边还是右边
			{
				min_cost_R = temp1;
				//NL += NP;
				temp_flag = 0;//重叠包围盒划分标志

			}
			else
			{
				min_cost_R = temp2;
				temp_flag = 1;
			}
			if (min_cost_R < min_cost)
			{
				min_cost = min_cost_R;
				Position->point[i] = (array + j)->bmax[i];//记录分割点位
				Position->PointsIndex =  NL ;//归属于左边的面元数量
				*np = NP;
				flag = i;//当在某维度上出现更小的代价时，记录其维度
				*box_flag = temp_flag;
			}
		}
	}
	return flag;//返回分割维度
}
/****************************************************************************************
函数功能：更新节点内的包围盒线索值等信息
参数意义：
返回值意义：选择的分割轴 0.1.2分别代表x.y.z轴
****************************************************************************************/
void Updata_Node(KD_Node *treenode, Prim_Box *array)
{
	for (int i = 0; i < 6; i++)//当前节点继承父节点的线索值,
	{
		treenode->LeftChild->rope[i] = treenode->rope[i];
		treenode->RightChild->rope[i] = treenode->rope[i];
	}
	//左右孩子继承父节点的包围盒的信息
	for (int i = 0; i < 3; i++)
	{
		treenode->LeftChild->box.bmin[i] = treenode->box.bmin[i];
		treenode->LeftChild->box.bmax[i] = treenode->box.bmax[i];
		treenode->RightChild->box.bmin[i] = treenode->box.bmin[i];
		treenode->RightChild->box.bmax[i] = treenode->box.bmax[i];
	}

	if (treenode->Split_Axis == -1)
	{
		printf("error when choose split");
		//StopBuild(treenode, array);
		return;
	}

	else if (treenode->Split_Axis == 0)
	{
		treenode->LeftChild->rope[1] = treenode->RightChild;
		treenode->RightChild->rope[0] = treenode->LeftChild;
		treenode->LeftChild->box.bmax[0] = treenode->SplitPos.point[treenode->Split_Axis];
		treenode->RightChild->box.bmin[0] = treenode->SplitPos.point[treenode->Split_Axis];
	}
	else if (treenode->Split_Axis == 1)
	{
		treenode->LeftChild->rope[3] = treenode->RightChild;
		treenode->RightChild->rope[2] = treenode->LeftChild;
		treenode->LeftChild->box.bmax[1] = treenode->SplitPos.point[treenode->Split_Axis];
		treenode->RightChild->box.bmin[1] = treenode->SplitPos.point[treenode->Split_Axis];
	}
	else if (treenode->Split_Axis == 2)
	{
		treenode->LeftChild->rope[5] = treenode->RightChild;
		treenode->RightChild->rope[4] = treenode->LeftChild;
		treenode->LeftChild->box.bmax[2] = treenode->SplitPos.point[treenode->Split_Axis];
		treenode->RightChild->box.bmin[2] = treenode->SplitPos.point[treenode->Split_Axis];
	}
	else printf("\nerror when get Split Axis\n");
}
//void Build_Tree(struct KD_Node * treenode,struct Prim_Box *array)
//{//建树过程 
//	int i;
//	DATA_Type max, min;
//	if (treenode->PrimCount  < Prim_Min || treenode->end - treenode->begin  < Prim_Min || treenode->Depth >= 20 || treenode->IsEmpty)
//	{
//		StopBuild(treenode, array);
//		return;
//	}
//	
//	else
//	{
//		treenode->PrimCount = treenode -> end - treenode -> begin + 1;//更新节点内三角面元数量
//		treenode->LeftChild = (KD_Node*)malloc(sizeof(KD_Node));
//		treenode->RightChild = (KD_Node*)malloc(sizeof(KD_Node));
//		treenode->LeftChild->IsLeaf = false;
//		treenode->RightChild->IsLeaf = false;
//		treenode->LeftChild->IsEmpty = false;
//		treenode->RightChild->IsEmpty = false;
//		treenode->LeftChild->Depth = treenode->Depth+1;//更新左孩子树高
//		treenode->RightChild->Depth = treenode->Depth+1;//更新右孩子树高
//		
//		for( i = 0; i < 6 ;i++)//当前节点继承父节点的线索值,
//		{
//			treenode->LeftChild->rope[i] = treenode->rope[i];
//			treenode->RightChild->rope[i] = treenode->rope[i];
//		}
////按照三个维度分别排序,选择分割轴并找出分割点
//
//
//		////执行的函数
//		treenode->Split_Axis = SAH_Choose_Split_Axis(treenode, array, treenode->begin, treenode->end, &treenode->SplitPos);
//		//Max_Min(array, treenode->begin, treenode->end, treenode->Split_Axis, &max, &min);
//
////左右孩子继承父节点的包围盒的信息
//		for (i = 0; i < 3; i++)
//		{
//			treenode->LeftChild->box.bmin[i] = treenode->box.bmin[i];
//			treenode->LeftChild->box.bmax[i] = treenode->box.bmax[i];
//			treenode->RightChild->box.bmin[i] = treenode->box.bmin[i];
//			treenode->RightChild->box.bmax[i] = treenode->box.bmax[i];
//		}
//		for (int i = 0; i < 3; i++)
//		{
//			if (treenode->box.bmax[i] < treenode->box.bmin[i])
//				printf("Box split error");
//		}
//		//QuickSort(array, treenode->begin, treenode->end, treenode->Split_Axis);//在退出前按照具有最小成本的维度重新排序
//
////判断按照选出的分割点进行分割，是否不可再分割
//		if (treenode->SplitPos.PointsIndex == 0)//分割位置在该节点的包围盒最左边面元
//		{
//			if (treenode->SplitPos.point[treenode->Split_Axis] == treenode->box.bmin[treenode->Split_Axis])//若最优分割位置不存在,停止分割
//			{
//				StopBuild(treenode, array);
//				return;
//			}
//			else //空节点j
//			{
//				treenode->LeftChild->IsEmpty = true;
//				treenode->LeftChild->PrimCount = 0;
//				treenode->RightChild->begin = treenode->begin;
//				treenode->RightChild->end = treenode->end;
//				treenode->RightChild->PrimCount = treenode->end - treenode->begin + 1;
//			}
//		}
//		else if (treenode->SplitPos.PointsIndex == treenode->PrimCount)
//		{
//			if (treenode->SplitPos.point[treenode->Split_Axis] == treenode->box.bmax[treenode->Split_Axis])
//			{
//					StopBuild(treenode, array);
//					return;
//			}
//			else
//			{
//				treenode->RightChild->IsEmpty = true;
//				treenode->RightChild->PrimCount = 0;
//				treenode->LeftChild->begin = treenode->begin;
//				treenode->LeftChild->end = treenode->end;
//				treenode->LeftChild->PrimCount = treenode->end - treenode->begin + 1;
//			}
//		}
//		else//更新子节点包含面元的起始范围
//		{
//			treenode->LeftChild->begin = treenode->begin;
//			treenode->LeftChild->end =treenode->begin + treenode->SplitPos.PointsIndex-1;
//			treenode->RightChild->begin = treenode->begin + treenode->SplitPos.PointsIndex;
//			treenode->RightChild->end = treenode->end;
//			treenode->LeftChild->PrimCount = treenode->LeftChild->end - treenode->LeftChild->begin+1;
//			treenode->RightChild->PrimCount = treenode->RightChild->end - treenode->RightChild->begin + 1;
//		}
//		/*==================
//		根据分割轴更新左右孩子的线索值和包围盒的值
//		===================*/
//		if (treenode->Split_Axis == -1)
//		{
//			StopBuild(treenode, array);
//			return;
//		}
//			//printf("\nerror when Choose Split Axis\n");
//		else if (treenode->Split_Axis == 0)
//		{
//			treenode->LeftChild->rope[1] = treenode->RightChild;
//			treenode->RightChild->rope[0] = treenode->LeftChild;
//			treenode->LeftChild->box.bmax[0] = treenode->SplitPos.point[treenode->Split_Axis];
//			treenode->RightChild->box.bmin[0] = treenode->SplitPos.point[treenode->Split_Axis];
//		}
//		else if (treenode->Split_Axis == 1)
//		{
//			treenode->LeftChild->rope[3] = treenode->RightChild;
//			treenode->RightChild->rope[2] = treenode->LeftChild;
//			treenode->LeftChild->box.bmax[1] = treenode->SplitPos.point[treenode->Split_Axis];
//			treenode->RightChild->box.bmin[1] = treenode->SplitPos.point[treenode->Split_Axis];
//		}
//		else if (treenode->Split_Axis == 2)
//		{
//			treenode->LeftChild->rope[5] = treenode->RightChild;
//			treenode->RightChild->rope[4] = treenode->LeftChild;
//			treenode->LeftChild->box.bmax[2] = treenode->SplitPos.point[treenode->Split_Axis];
//			treenode->RightChild->box.bmin[2] = treenode->SplitPos.point[treenode->Split_Axis];
//		}
//		else printf("\nerror when get Split Axis\n");
//
//	}
//
//	Build_Tree(treenode->LeftChild, array);//递归出左右孩子
//	Build_Tree(treenode->RightChild, array);
//}
void Build_BigNode(struct KD_Node * treenode, struct Prim_Box *array ,struct Prim_Box *out_array,int arr_length,int *out_arr_length)
{
	int i,nl,nr,np,box_flag;
	DATA_Type max, min;
	Prim_Box *arr_left, *arr_right;
	arr_left = NULL;
	arr_right = NULL;
	if (treenode == NULL)
	{
		free(treenode);
		return;
	}
	if ((arr_length  < Prim_Min) || (treenode->IsEmpty) || (treenode->Depth == Max_Depth))
	{
		StopBuild(treenode, array,out_array,arr_length,out_arr_length);
		return;
	}
	else
	{
		//treenode->PrimCount = arr_length;//更新节点内三角面元数量
		treenode->LeftChild = (KD_Node*)malloc(sizeof(KD_Node));
		treenode->RightChild = (KD_Node*)malloc(sizeof(KD_Node));
		treenode->LeftChild->IsLeaf = false;
		treenode->RightChild->IsLeaf = false;
		treenode->LeftChild->IsEmpty = false;
		treenode->RightChild->IsEmpty = false;
		treenode->LeftChild->Depth = treenode->Depth + 1;//更新左孩子树高
		treenode->RightChild->Depth = treenode->Depth + 1;//更新右孩子树高
		//中分法选分割轴
		if (treenode->PrimCount > BigNodeNum)
		{
			if (!Cut_Blank(treenode, &array,&arr_left,&arr_right,arr_length,&nl,&nr))
			{
				treenode->Split_Axis = Mid_Choose_Split_Axis(treenode, array, arr_length, &treenode->SplitPos);
				
				Count_Prim_Mid(array, arr_length, treenode->SplitPos.point[treenode->Split_Axis], &nl, &nr, treenode->Split_Axis);
				
				arr_left = (Prim_Box*)malloc(nl*sizeof(Prim_Box));
				arr_right = (Prim_Box*)malloc(nr*sizeof(Prim_Box));

				Resize_Array_Mid(array, arr_left, arr_right, arr_length, nl, nr, treenode->SplitPos.point[treenode->Split_Axis], treenode->Split_Axis);
				treenode->LeftChild->PrimCount = nl;
				treenode->RightChild->PrimCount = nr;
				Updata_Node(treenode, array);
				//free(array);
				//array=NULL;
			}
		}
		else
		{
			treenode->Split_Axis = SAH_Choose_Split_Axis(treenode, array, arr_length, &box_flag,&treenode->SplitPos,&np);
			//判断按照选出的分割点进行分割，是否不可再分割

			//if ((treenode->SplitPos.PointsIndex + np) == 0)//分割位置在该节点的包围盒最左边面元
			//{
				if (treenode->SplitPos.point[treenode->Split_Axis] == treenode->box.bmin[treenode->Split_Axis])//若最优分割位置不存在,停止分割
				{
					StopBuild(treenode, array,out_array,arr_length,out_arr_length);
					free(treenode->LeftChild);
					free(treenode->RightChild);
					treenode->LeftChild = NULL;
					treenode->RightChild = NULL;
					return;
				}
			//	else //空节点j
			//	{
			//		treenode->LeftChild->IsEmpty = true;
			//		treenode->LeftChild->PrimCount = 0;
			//		treenode->RightChild->PrimCount = arr_length;
			//		nl = 0;
			//		nr = arr_length;
			//		arr_right = array;
			//	}
			//}
			//else if (((treenode->SplitPos.PointsIndex + np) == treenode->PrimCount)) 
			//{
				else if (treenode->SplitPos.point[treenode->Split_Axis] == treenode->box.bmax[treenode->Split_Axis])
				{
					StopBuild(treenode, array, out_array, arr_length, out_arr_length);
					free(treenode->LeftChild);
					free(treenode->RightChild);
					treenode->LeftChild = NULL;
					treenode->RightChild = NULL;
					return;
				}
				//else
				//{
				//	treenode->RightChild->IsEmpty = true;
				//	treenode->RightChild->PrimCount = 0;
				//	treenode->LeftChild->PrimCount = arr_length;
				//	nl = arr_length;
				//	nr = 0;
				//	arr_left = array;
				//}
			//}
			else//更新子节点包含面元的起始范围
			{
				//treenode->LeftChild->begin = treenode->begin;
				//treenode->LeftChild->end = treenode->begin + treenode->SplitPos.PointsIndex - 1;
				//treenode->RightChild->begin = treenode->begin + treenode->SplitPos.PointsIndex;
				//treenode->RightChild->end = treenode->end;
				//Count_Prim(array, arr_length, treenode->SplitPos.point[treenode->Split_Axis], &nl, &nr, &np, &treenode->Split_Axis);
				nl = treenode->SplitPos.PointsIndex;
				nr = arr_length - nl - np;
				arr_left = (Prim_Box*)malloc((nl + np)*sizeof(Prim_Box));
				arr_right = (Prim_Box*)malloc((nr + np)*sizeof(Prim_Box));
				Resize_Array(array, arr_left, arr_right, arr_length, (nl+np), (nr+np), treenode->SplitPos.point[treenode->Split_Axis], treenode->Split_Axis);
				treenode->LeftChild->PrimCount = (nl + np);
				treenode->RightChild->PrimCount = (nr + np);
				
			}
			/*==================
			根据分割轴更新左右孩子的线索值和包围盒的值
			===================*/
			Updata_Node(treenode, array);
			//free(array);
			//array=NULL;
		}
	}
	Build_BigNode(treenode->LeftChild, arr_left, out_array, treenode->LeftChild->PrimCount, out_arr_length);//递归出左右孩子
	Build_BigNode(treenode->RightChild, arr_right, out_array, treenode->RightChild->PrimCount, out_arr_length);
	//free(array);
	//free(arr_right);
	//array= NULL;
	//arr_right = NULL;
}

/***************************************************
函数功能：定义叶节点，算出节点内的三角面元的数量，
并列出节点内三角面元的列表将属于该节点
的三角面元的值赋给节点内部的结构体并进行保存
参数意义：
****************************************************/
void StopBuild(KD_Node * treenode, Prim_Box *array, Prim_Box *out_array,int arr_length,int *out_arr_length)
{
	if (treenode->IsEmpty == true)
	{
		//free(treenode->LeftChild);
		//free(treenode->RightChild);
		treenode->LeftChild = NULL;
		treenode->RightChild = NULL;
		treenode->PrimCount = 0;
		return;
	}
	else
	{
		treenode->IsLeaf = true;
		treenode->PrimCount = arr_length;
		for (int i = 0; i < arr_length; i++)
		{
			*(out_array + i + *out_arr_length) = *(array + i);
		}
		treenode->begin = *out_arr_length;
		treenode->end = *out_arr_length + arr_length - 1;
		*out_arr_length = *out_arr_length + arr_length;
		free(array);
		//free(treenode->LeftChild);
		//free(treenode->RightChild);
		treenode->LeftChild = NULL;
		treenode->RightChild = NULL;
		array = NULL;
	
	}
	return;
		
}
void Optimization_Rope(KD_Node *treenode)
{
	if (treenode->LeftChild == NULL || treenode->LeftChild->IsEmpty == true || treenode->LeftChild->IsLeaf == true 
		|| treenode->RightChild == NULL || treenode->RightChild->IsEmpty == true || treenode->RightChild->IsLeaf == true)
		return;
	if (treenode->LeftChild->Split_Axis == treenode->RightChild->Split_Axis && treenode->Split_Axis == treenode->LeftChild->Split_Axis)//父节点左右子节点分割轴相同
	{
		treenode->LeftChild->RightChild->rope[2 * treenode->Split_Axis + 1] = treenode->RightChild->LeftChild;
		treenode->RightChild->LeftChild->rope[2 * treenode->Split_Axis] = treenode->LeftChild->RightChild;
	}
	else if (treenode->LeftChild->Split_Axis == treenode->RightChild->Split_Axis && treenode->Split_Axis != treenode->LeftChild->Split_Axis)
	{
		if (treenode->RightChild->SplitPos.point[treenode->RightChild->Split_Axis] < treenode->LeftChild->SplitPos.point[treenode->LeftChild->Split_Axis])//右子节点分隔位置小于左子节点分割位置
		{
			treenode->LeftChild->RightChild->rope[2 * treenode->Split_Axis + 1] = treenode->RightChild->RightChild;
			treenode->RightChild->LeftChild->rope[2 * treenode->Split_Axis] = treenode->LeftChild->LeftChild;
		}
		else if (treenode->RightChild->SplitPos.point[treenode->RightChild->Split_Axis] > treenode->LeftChild->SplitPos.point[treenode->LeftChild->Split_Axis])
		{
			treenode->LeftChild->LeftChild->rope[2 * treenode->Split_Axis + 1] = treenode->RightChild->LeftChild;
			treenode->RightChild->RightChild->rope[2 * treenode->Split_Axis] = treenode->LeftChild->RightChild;
		}
	}
	else if (treenode->Split_Axis == treenode->LeftChild->Split_Axis && treenode->LeftChild->Split_Axis != treenode->RightChild->Split_Axis)
	{
		treenode->RightChild->LeftChild->rope[2 * treenode->Split_Axis] = treenode->LeftChild->RightChild;
		treenode->RightChild->RightChild->rope[2 * treenode->Split_Axis] = treenode->LeftChild->RightChild;
	}
	else if (treenode->Split_Axis == treenode->RightChild->Split_Axis && treenode->LeftChild->Split_Axis != treenode->RightChild->Split_Axis)
	{
		treenode->LeftChild->RightChild->rope[2 * treenode->Split_Axis + 1] = treenode->RightChild->LeftChild;
		treenode->LeftChild->LeftChild->rope[2 * treenode->Split_Axis + 1] = treenode->RightChild->LeftChild;
	}
	Optimization_Rope(treenode->LeftChild);
	Optimization_Rope(treenode->RightChild);
}
//
//void Destroy_Tree(KD_Node *root)
//{
//	if(root)
//	{
//		Destroy_Tree(root->LeftChild);
//		Destroy_Tree(root->RightChild);
//	}
//	free(root);
//	root=NULL;
//}




