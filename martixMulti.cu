
//矩阵相乘cpu版本
//输入：x,y相乘的两矩阵，一维数组。	rowNum，colNum分别为两矩阵的维数。 result为矩阵相乘结果，一维数组
//返回值：-1失败   1成功
int martixMulti(float* x,int rowNum_X,int colNum_X,float* y,int rowNum_Y,int colNum_Y,float* result)
{
	if(colNum_X != rowNum_Y)
		return -1;
	
	memset(result,0,sizeof(float) * rowNum_X*colNum_Y);

	for(int i=0;i<rowNum_X;i++)
	{
		for(int j=0;j<colNum_Y;j++)
		{
			for(int idx=0;idx<colNum_X;idx++)
			{
				result[i*colNum_Y +j] += x[i*colNum_X + idx] * y[idx*colNum_Y + j];
			}
		}
	}
	return 1;
}
