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
�������ܣ���������Ԫ����Ԥ����ÿ����Ԫ��ʹ�ð�Χ�а�����������SAH����
�������壺1��Ҫ�򿪵��ļ�ָ��
2��������������Ԫ������
3�����涥����Ϣ
4�����潨������Ԫ��Χ��
********************************************************************************************************/
void Preprocessing_Triangles(FILE *fp, Triangle *Triangles, Element *Points, int *Node_Num,int *Triangle_Num)
{
	int i, j;
	int temp;//�������λ��������

/*-----��ȡdat�ļ�------*/
	int node_num;
	int face_num;
	int index = 0,num = 0;
	int tempnum1,tempnum2,tempnum3;

	fscanf(fp,"%d",&num);
	printf("node_num=%d\n",num);						// ���ڵ���;
	float *pointx = new float[num];					// �ڵ��x����
	float *pointy = new float[num];					// �ڵ��y����
	float *pointz = new float[num];					// �ڵ��z����
	while (num)
	{

	        fscanf(fp,"%d",&index);

		///////�������������/////
		fscanf(fp,"%f",&(pointx[index-1]));			// ��x����;     
		fscanf(fp,"%f",&(pointy[index-1]));			// ��y����;
		fscanf(fp,"%f",&(pointz[index-1]));			// ��z����;
                //////////////////////////
		
		///////��(benchmark)����/////
		//fscanf(fp,"%d",&index);
		//fscanf(fp,"%f",&(pointz[index-1]));			// ��x����;     
		//fscanf(fp,"%f",&(pointy[index-1]));			// ��y����;
		//fscanf(fp,"%f",&(pointx[index-1]));			// ��z����;
		//pointx[index-1] = pointx[index-1]  + 24.0*1.5;
		//pointx[index-1] = pointx[index-1] ;
		//pointy[index-1] = pointy[index-1] ;
		//pointz[index-1] = pointz[index-1] ;
		//////////////////////////////

		num--;
	}
	node_num = index;                                               //�ڵ����

	fscanf(fp,"%d",&num);	
	printf("ele_num=%d\n",num);					// ����Ԫ��;
	int *faceOne = new int[num];					// �ڵ��x����
	int *faceTwo = new int[num];					// �ڵ��y����
	int *faceThree = new int[num];					// �ڵ��z����
	while (num)
	{
		fscanf(fp,"%d",&index);
		fscanf(fp,"%d",&(faceOne[index-1]));			// ��x����;     
		fscanf(fp,"%d",&(faceTwo[index-1]));			// ��y����;
		fscanf(fp,"%d",&(faceThree[index-1]));			// ��z����;
		num--;
	}
	
	face_num = index;                                                //��Ԫ����




	for (i = 0; i < node_num; i++)
	{
		
			/*(Points + i)->point[0] = pointx[i];
			(Points + i)->point[1] = pointy[i];
			(Points + i)->point[2] = pointz[i];*/

	        Points[i].point[0] = pointx[i];
			Points[i].point[1] = pointy[i];
			Points[i].point[2] = pointz[i];

		//fscanf(fp,"%d",&(Points + i)->PointsIndex);//��������б��
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
	int temp;//�������λ��������


	for (i = 0; i < Points_Num; i++)
	{
		for (j = 0; j < 3; j++)
			fscanf(fp, "%f", &(Points + i)->point[j]);
		//fscanf(fp,"%d",&(Points + i)->PointsIndex);//��������б��
		(Points + i)->PointsIndex = i;
		//printf("\nRead  Point%d", (Points + i)->PointsIndex);
	}

	for (i = 0; i<Triangle_Num; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(fp,"%x",&(Triangles + i)->Points[j]);//������Ԫ�ڴ洢�Ķ�����
			(Triangles + i)->Points[j] = (Triangles + i)->Points[j] -1 ;//���ڶ�������ֵ��1��ʼ�����Ҫ��һ
		}
		fscanf(fp, "%d", &temp);//������Ԫ�ڴ洢�Ķ�����
		fscanf(fp, "%d", &temp);//������Ԫ�ڴ洢�Ķ�����
		(Triangles + i)->TriangleIndex = i;
		//fscanf(fp,"%d",&(Triangles + i)->TriangleIndex);//�����α��
		//printf("\nRead  Triangle%d", (Triangles + i)->TriangleIndex);
	}

	Sort_Box(Points, Triangles, Prim, Triangle_Num);
}*/
/*******************************************************
�������ܣ��������������Ԫ�İ�Χ�е��������������ֵ
��������:1����������������Ԫ������
2���������ж��������
3����������Ԫ��Χ�й��ɵ�����
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
			//else//�������δ�ֱ��������ʱ������һ�����
			//{
			//	(prim + i)->bmax[j] = Max+Lambda/2;
			//	(prim + i)->bmin[j] = Min - Lambda/2;
			//}
		}
		(prim + i)->Box_Index = i;
		//ÿ����Χ���������ε�����
		(prim + i)->bmid[0] = ((points + (triangle + i)->Points[0])->point[0] + (points + (triangle + i)->Points[1])->point[0] + (points + (triangle + i)->Points[2])->point[0]) / 3;
		(prim + i)->bmid[1] = ((points + (triangle + i)->Points[0])->point[1] + (points + (triangle + i)->Points[1])->point[1] + (points + (triangle + i)->Points[2])->point[1]) / 3;
		(prim + i)->bmid[2] = ((points + (triangle + i)->Points[0])->point[2] + (points + (triangle + i)->Points[1])->point[2] + (points + (triangle + i)->Points[2])->point[2]) / 3;
	}

}
/****************************************************************************************
�������ܣ��Խڵ���г�ʼ��
�������壺1����ʼ�ڵ�
2��������Ԫ��Χ�й��ɵ�����
����ֵ���壺ѡ��ķָ��� 0.1.2�ֱ����x.y.z��
****************************************************************************************/
void KD_Node_init(struct KD_Node *kd_node, Prim_Box *array,int arr_length)
{//�ڵ��ʼ��
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
	//����ʼ�ڵ�İ�Χ��ȷ����Χ�����ҵ���Χ������ԭ���������Զ�ĵ�
	Max_Min(array, arr_length, 0, &kd_node->box.bmax[0], &kd_node->box.bmin[0]);
	Max_Min(array, arr_length, 1, &kd_node->box.bmax[1], &kd_node->box.bmin[1]);
	Max_Min(array, arr_length, 2, &kd_node->box.bmax[2], &kd_node->box.bmin[2]);
	return;
}
/******************************** 
 *��������swap 
 *���ã����������ṹ���ֵ 
 *�����������������ṹ�� 
 *����ֵ���� 
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
 *��������quicksort 
 *���ã����������㷨����С��������
 *������ ������Ľṹ������ �������е���ʼλ�ã������н���λ�ã���Ҫ�����ά��
 *����ֵ���� 
 ************************************/  
void QuickSort(Prim_Box* array,  int begin, int end,int dim) 
{  
    int i, j;  
    if(begin < end)  
    {  
        i = begin + 1;  // ��array[begin]��Ϊ��׼������˴�array[begin+1]��ʼ���׼���Ƚϣ�  
        j = end;        // array[end]����������һλ  
            
        while(i < j)  
        {  
            if(  (array+i)->bmid[dim] > (array+begin) ->bmid[dim])  // ����Ƚϵ�����Ԫ�ش��ڻ�׼�����򽻻�λ�á�  
            {  
                swap(&array[i],&array[j]);  // ���������ṹ��
                j--;  
            }  
            else  
            {  
                i++;  // �����������һλ���������׼���Ƚϡ�  
            }  
        }  
		/* ����whileѭ����i = j�� 
         * ��ʱ���鱻�ָ����������  -->  array[begin+1] ~ array[i-1] < array[begin] 
         *                           -->  array[i+1] ~ array[end] > array[begin] 
         * ���ʱ������array�ֳ��������֣��ٽ�array[i]��array[begin]���бȽϣ�����array[i]��λ�á� 
         * ���array[i]��array[begin]���������������ָ�ֵ������Դ����ƣ�ֱ�����i = j�������������˳��� 
         */  
        if( (array+i)->bmid[dim] >= (array+begin)->bmid[dim])  // �������Ҫȡ�ȡ�>=������������Ԫ������ͬ��ֵʱ������ִ���  
        {  
            i--;  
        }  
        swap(&array[begin], &array[i]);  // ����array[i]��array[begin]  
        QuickSort(array, begin, i,dim);  
        QuickSort(array, j, end,dim);  
    }  
}  

/****************************************************************************************
�������ܣ�ʹ���зַ�ѡ��ָ���
�������壺1���ɰ�Χ�й��ɵ�����
2��3��begin��end��ʾ���������������е���ʼ��
4��mid��ʾ�ָ���ֵ
����ֵ���壺ѡ��ķָ��� 0.1.2�ֱ����x.y.z��
****************************************************************************************/
int Mid_Choose_Split_Axis(KD_Node *treenode, Prim_Box *array, int arr_length, Element* mid)
{
	//���������ϵĳ���
	DATA_Type x_dim = 0;
	DATA_Type y_dim = 0;
	DATA_Type z_dim = 0;
	//���������ϵ��м�ֵ�����ָ��
	DATA_Type x_mid = 0;
	DATA_Type y_mid = 0;
	DATA_Type z_mid = 0;
	int nl,nr;
	//��������������������
	//int x_y_z = 0;//0��ʾx����1��ʾy��2��ʾz
	/***********************************************************
	������ά���ϵĳ��ȣ����ҳ�����Ǹ����򣬲���¼�е�ֵ
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
�������ܣ��޳���Χ���еĿհײ���

*********************************************************/
bool Cut_Blank(KD_Node *treenode, Prim_Box **array, Prim_Box **arr_left, Prim_Box **arr_right, int arr_length, int *nl, int *nr)
{
	DATA_Type max, min, blank, length,empty_rate;
	empty_rate = 0.20;
	for (int i = 0; i < 3; i++)//����ά�����ҿհײ��֣���������ֵʱ���зָ�
	{
		Max_Min(*array, arr_length, i, &max, &min);
		length = max - min;
		blank = fabs(max - treenode->box.bmax[i]);//x�����Һ���Ϊ�սڵ�
		if (blank / length > empty_rate)
		{
			treenode->RightChild->IsEmpty = true;
			treenode->RightChild->PrimCount = 0;

			treenode->LeftChild->PrimCount = arr_length;
			//arr_left = (Prim_Box*)malloc(arr_length*sizeof(Prim_Box));
			*arr_left = *array;//�Һ���Ϊ�գ�����ǰ����ֱ��ת������
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
�������ܣ��ҳ������е������Сֵ
�������壺
1��������
2�����������
3�������ҷָ���
4�����ֵ
5����Сֵ
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
�������ܣ�����ָ�ƽ�������Լ����ָ�ƽ�洩�����Ӱ�Χ�е�����
�������壺1���ɰ�Χ�й��ɵ�����
2��3��begin��end��ʾ���������������е���ʼ��
4��Poisition�ڰ����˷ָ��λ���Լ��ָ�ƽ�����ڰ�Χ�еı��
5���ýڵ��а�Χ�����ֵС�ڷָ�����Ԫ����
6���ýڵ��а�Χ�����ֵ���ڷָ�����Ԫ����
7���ýڵ��а�Χ�б��ָ�ƽ��ָ����Ԫ����
8���ָ���
����ֵ���壺�޷���ֵ
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
			(*nl) += ((itr->bmax[*Axis]) <= Split_Pos ? 1 : 0);//���޸�
			(*nr) += ((itr->bmin[*Axis]) >= Split_Pos ? 1 : 0);
		}

	}
	*np = arr_length - *nr - *nl ;
	return;
}
/****************************************************************************************
�������ܣ�����ָ�ƽ�������Լ����ָ�ƽ�洩�����Ӱ�Χ�е�����
�������壺1���ɰ�Χ�й��ɵ�����
2��3��begin��end��ʾ���������������е���ʼ��
4��Poisition�ڰ����˷ָ��λ���Լ��ָ�ƽ�����ڰ�Χ�еı��
5���ýڵ��а�Χ�����ֵС�ڷָ�����Ԫ����
6���ýڵ��а�Χ�����ֵ���ڷָ�����Ԫ����
7���ýڵ��а�Χ�б��ָ�ƽ��ָ����Ԫ����
8���ָ���
����ֵ���壺�޷���ֵ
****************************************************************************************/
void Count_Prim_Mid(Prim_Box *array, int arr_length, DATA_Type Split_Pos, int *nl, int *nr, int Axis)
{
	*nr = 0;
	*nl = 0;
	for (Prim_Box *itr = array ; itr != (array + arr_length); itr++)
	{
		(*nl) += ((itr->bmin[Axis]) <= Split_Pos ? 1 : 0);//���޸�
		(*nr) += ((itr->bmax[Axis]) >= Split_Pos ? 1 : 0);
	}
	return;
}
/****************************************************************************************
�������ܣ������鰴���ʷ�ƽ�������з֣�������ԭ����λ��
�������壺
1��������
2��3�����������Ԫ�������е���ʼλ������ֹλ��
4��5���ֱ��������ҽڵ���Ԫ������
6���ص���Ԫ������߻����ұߵı�־λ 0����� 1���ұ�
7���ָ�λ��
8���ָ���
****************************************************************************************/
void Resize_Array(Prim_Box *array, Prim_Box *arr_left, Prim_Box *arr_right, int arr_length, int nl, int nr, int flag, DATA_Type Split_Pos, int Axis)
{
	//struct Prim_Box *left = (Prim_Box*)malloc(nl*(sizeof(Prim_Box)));
	//struct Prim_Box *right = (Prim_Box*)malloc(nr*(sizeof(Prim_Box)));
	//memset(left, 0, nl*sizeof(struct Prim_Box));//��ʼ��Ϊ0
	//memset(right, 0, nr*sizeof(struct Prim_Box));//��ʼ��Ϊ0

	if (!flag)//�ص���Ԫ�����
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
	//memset(left, 0, nl*sizeof(struct Prim_Box));//��ʼ��Ϊ0
	//memset(right, 0, nr*sizeof(struct Prim_Box));//��ʼ��Ϊ0
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
	//memset(left, 0, nl*sizeof(struct Prim_Box));//��ʼ��Ϊ0
	//memset(right, 0, nr*sizeof(struct Prim_Box));//��ʼ��Ϊ0
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
�������ܣ�ʹ��SAH�㷨ѡ��ָ���
�������壺1���ɰ�Χ�й��ɵ�����
2��3��begin��end��ʾ���������������е���ʼ��
4��Poisition�ڰ����˷ָ��λ���Լ��ָ�ƽ�����ڰ�Χ�еı��
����ֵ���壺ѡ��ķָ��� 0.1.2�ֱ����x.y.z��
****************************************************************************************/
int SAH_Choose_Split_Axis(KD_Node *treenode, Prim_Box *array, int arr_length, int *box_flag,struct Element *Position,int *np)
{
	int flag = -1;//ѡȡ��������
	//int box_flag = 0;
	int temp_flag = 0;
	DATA_Type temp1, temp2;//ѡ��ڵ����Ҳ�ı�־
	int c_hit, c_walk,NL,NR,NP;
	c_hit = 19;
	c_walk = 9;//??????????????????????????????????????????
	//int alpha_k = 5;//alpha_k=c_hit/c_walk;
	DATA_Type cost = 0;
	DATA_Type min_cost, min_cost_L, min_cost_R;//ѡ����ʷ�����ѡ�а�Χ�е���߻����ұ�
	DATA_Type Surface,Surface_L,Surface_R;
	DATA_Type H, L, W,
					   H_L,L_L,W_L,
					   H_R,L_R,W_R;
	DATA_Type max, min;
		//��ǰ�ڵ��Χ�еĳ����
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
	min_cost = DATA_Type(c_hit * 10 * (arr_length)+c_walk);//��ʼ����С�Ĵ���Ϊ���ֵ
	for (int i = 0; i < 3; i++)
	{
		min_cost_L = min_cost;
		min_cost_R = min_cost;
		//QuickSort(array, begin, end, i);
		for (int j = 0; j <arr_length; j++)
		{
			if ((array + j)->bmin[i] < treenode->box.bmin[i] || (array + j)->bmax[i] > treenode->box.bmax[i])
				continue;
			switch (i)//�԰�Χ�е���ƽ����Ϊ�ʷ��棬��ѡ���ʷֽڵ�İ�Χ�б������ҽڵ���
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
			if (temp1 <= temp2 )//�ص��İ�Χ���ǻ�����߻����ұ�
			{
				min_cost_L = temp1;
				//NL += NP;
				temp_flag = 0;//�ص���Χ�л��ֱ�־
			}
			else
			{
				min_cost_L = temp2;
				temp_flag = 1;
			}
			if (min_cost_L < min_cost)
			{
				Position->point[i] = (array + j)->bmin[i];//��¼�ָ��λ
				Position->PointsIndex = NL ;//���ص��ڵ���ൽ��ߣ����Χ�д�С��ԭ����һ�£��߽���
				*np = NP;
				min_cost = min_cost_L;
				flag = i;//����ĳά���ϳ��ָ�С�Ĵ���ʱ����¼��ά��
				*box_flag = temp_flag;
			}
			//min_cost_L = min_cost_L <= cost ? min_cost_L : cost;
			switch (i)//�԰�Χ�е���ƽ����Ϊ�ʷ��棬��ѡ���ʷֽڵ�İ�Χ�б�������ڵ���
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
			if (temp1 <= temp2)//�ص��İ�Χ���ǻ�����߻����ұ�
			{
				min_cost_R = temp1;
				//NL += NP;
				temp_flag = 0;//�ص���Χ�л��ֱ�־

			}
			else
			{
				min_cost_R = temp2;
				temp_flag = 1;
			}
			if (min_cost_R < min_cost)
			{
				min_cost = min_cost_R;
				Position->point[i] = (array + j)->bmax[i];//��¼�ָ��λ
				Position->PointsIndex =  NL ;//��������ߵ���Ԫ����
				*np = NP;
				flag = i;//����ĳά���ϳ��ָ�С�Ĵ���ʱ����¼��ά��
				*box_flag = temp_flag;
			}
		}
	}
	return flag;//���طָ�ά��
}
/****************************************************************************************
�������ܣ����½ڵ��ڵİ�Χ������ֵ����Ϣ
�������壺
����ֵ���壺ѡ��ķָ��� 0.1.2�ֱ����x.y.z��
****************************************************************************************/
void Updata_Node(KD_Node *treenode, Prim_Box *array)
{
	for (int i = 0; i < 6; i++)//��ǰ�ڵ�̳и��ڵ������ֵ,
	{
		treenode->LeftChild->rope[i] = treenode->rope[i];
		treenode->RightChild->rope[i] = treenode->rope[i];
	}
	//���Һ��Ӽ̳и��ڵ�İ�Χ�е���Ϣ
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
//{//�������� 
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
//		treenode->PrimCount = treenode -> end - treenode -> begin + 1;//���½ڵ���������Ԫ����
//		treenode->LeftChild = (KD_Node*)malloc(sizeof(KD_Node));
//		treenode->RightChild = (KD_Node*)malloc(sizeof(KD_Node));
//		treenode->LeftChild->IsLeaf = false;
//		treenode->RightChild->IsLeaf = false;
//		treenode->LeftChild->IsEmpty = false;
//		treenode->RightChild->IsEmpty = false;
//		treenode->LeftChild->Depth = treenode->Depth+1;//������������
//		treenode->RightChild->Depth = treenode->Depth+1;//�����Һ�������
//		
//		for( i = 0; i < 6 ;i++)//��ǰ�ڵ�̳и��ڵ������ֵ,
//		{
//			treenode->LeftChild->rope[i] = treenode->rope[i];
//			treenode->RightChild->rope[i] = treenode->rope[i];
//		}
////��������ά�ȷֱ�����,ѡ��ָ��Ტ�ҳ��ָ��
//
//
//		////ִ�еĺ���
//		treenode->Split_Axis = SAH_Choose_Split_Axis(treenode, array, treenode->begin, treenode->end, &treenode->SplitPos);
//		//Max_Min(array, treenode->begin, treenode->end, treenode->Split_Axis, &max, &min);
//
////���Һ��Ӽ̳и��ڵ�İ�Χ�е���Ϣ
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
//		//QuickSort(array, treenode->begin, treenode->end, treenode->Split_Axis);//���˳�ǰ���վ�����С�ɱ���ά����������
//
////�жϰ���ѡ���ķָ����зָ�Ƿ񲻿��ٷָ�
//		if (treenode->SplitPos.PointsIndex == 0)//�ָ�λ���ڸýڵ�İ�Χ���������Ԫ
//		{
//			if (treenode->SplitPos.point[treenode->Split_Axis] == treenode->box.bmin[treenode->Split_Axis])//�����ŷָ�λ�ò�����,ֹͣ�ָ�
//			{
//				StopBuild(treenode, array);
//				return;
//			}
//			else //�սڵ�j
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
//		else//�����ӽڵ������Ԫ����ʼ��Χ
//		{
//			treenode->LeftChild->begin = treenode->begin;
//			treenode->LeftChild->end =treenode->begin + treenode->SplitPos.PointsIndex-1;
//			treenode->RightChild->begin = treenode->begin + treenode->SplitPos.PointsIndex;
//			treenode->RightChild->end = treenode->end;
//			treenode->LeftChild->PrimCount = treenode->LeftChild->end - treenode->LeftChild->begin+1;
//			treenode->RightChild->PrimCount = treenode->RightChild->end - treenode->RightChild->begin + 1;
//		}
//		/*==================
//		���ݷָ���������Һ��ӵ�����ֵ�Ͱ�Χ�е�ֵ
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
//	Build_Tree(treenode->LeftChild, array);//�ݹ�����Һ���
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
		//treenode->PrimCount = arr_length;//���½ڵ���������Ԫ����
		treenode->LeftChild = (KD_Node*)malloc(sizeof(KD_Node));
		treenode->RightChild = (KD_Node*)malloc(sizeof(KD_Node));
		treenode->LeftChild->IsLeaf = false;
		treenode->RightChild->IsLeaf = false;
		treenode->LeftChild->IsEmpty = false;
		treenode->RightChild->IsEmpty = false;
		treenode->LeftChild->Depth = treenode->Depth + 1;//������������
		treenode->RightChild->Depth = treenode->Depth + 1;//�����Һ�������
		//�зַ�ѡ�ָ���
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
			//�жϰ���ѡ���ķָ����зָ�Ƿ񲻿��ٷָ�

			//if ((treenode->SplitPos.PointsIndex + np) == 0)//�ָ�λ���ڸýڵ�İ�Χ���������Ԫ
			//{
				if (treenode->SplitPos.point[treenode->Split_Axis] == treenode->box.bmin[treenode->Split_Axis])//�����ŷָ�λ�ò�����,ֹͣ�ָ�
				{
					StopBuild(treenode, array,out_array,arr_length,out_arr_length);
					free(treenode->LeftChild);
					free(treenode->RightChild);
					treenode->LeftChild = NULL;
					treenode->RightChild = NULL;
					return;
				}
			//	else //�սڵ�j
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
			else//�����ӽڵ������Ԫ����ʼ��Χ
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
			���ݷָ���������Һ��ӵ�����ֵ�Ͱ�Χ�е�ֵ
			===================*/
			Updata_Node(treenode, array);
			//free(array);
			//array=NULL;
		}
	}
	Build_BigNode(treenode->LeftChild, arr_left, out_array, treenode->LeftChild->PrimCount, out_arr_length);//�ݹ�����Һ���
	Build_BigNode(treenode->RightChild, arr_right, out_array, treenode->RightChild->PrimCount, out_arr_length);
	//free(array);
	//free(arr_right);
	//array= NULL;
	//arr_right = NULL;
}

/***************************************************
�������ܣ�����Ҷ�ڵ㣬����ڵ��ڵ�������Ԫ��������
���г��ڵ���������Ԫ���б����ڸýڵ�
��������Ԫ��ֵ�����ڵ��ڲ��Ľṹ�岢���б���
�������壺
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
	if (treenode->LeftChild->Split_Axis == treenode->RightChild->Split_Axis && treenode->Split_Axis == treenode->LeftChild->Split_Axis)//���ڵ������ӽڵ�ָ�����ͬ
	{
		treenode->LeftChild->RightChild->rope[2 * treenode->Split_Axis + 1] = treenode->RightChild->LeftChild;
		treenode->RightChild->LeftChild->rope[2 * treenode->Split_Axis] = treenode->LeftChild->RightChild;
	}
	else if (treenode->LeftChild->Split_Axis == treenode->RightChild->Split_Axis && treenode->Split_Axis != treenode->LeftChild->Split_Axis)
	{
		if (treenode->RightChild->SplitPos.point[treenode->RightChild->Split_Axis] < treenode->LeftChild->SplitPos.point[treenode->LeftChild->Split_Axis])//���ӽڵ�ָ�λ��С�����ӽڵ�ָ�λ��
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




