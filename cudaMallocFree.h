

//void MallocOnGPU_getWidthHeight(int st, float fi, int Node_Num, Element* points, Radius** d_direction, Element** d_point, Radius** d_spherical, float** d_st_min, float** d_fi_max, float** d_st_max, float** d_fi_min);

void MallocOnGPU(int dwidth, int dheight, Direction** d_rays1, Square** d_squares1, Direction** d_rays2, Square** d_squares2,
	RayBeamInfo** d_effrays, Vector** d_center, Vector** d_axis, MatStruct** d_transMat, ReimOutput** d_reim,
	float** d_sum_re, float** d_sum_im, float** d_sum_sre,float** d_sum_sim, int** d_DivRayTubeNum, int** d_sum_gmem, int** d_sum_Gmem, int** d_squares_pred);

void MemsetOnGPU3(int dwidth, int dheight, RayBeamInfo** d_effrays, Vector** d_center, Vector** d_axis, MatStruct** d_transMat, ReimOutput** d_reim, float** d_sum_re, float** d_sum_im);

void MemsetOnGPU1(int dwidth, int dheight, Direction** d_rays1, Square** d_squares1, RayBeamInfo** d_effrays, Vector** d_center, Vector** d_axis,
	MatStruct** d_transMat, ReimOutput** d_reim, float** d_sum_re, float** d_sum_im,
	int** d_DivRayTubeNum, int** d_sum_gmem, int** d_sum_Gmem, int** d_squares_pred);

void MemsetOnGPU2(int dwidth, int dheight, Direction** d_rays2, Square** d_squares1, RayBeamInfo** d_effrays, Vector** d_center, Vector** d_axis,
	MatStruct** d_transMat, ReimOutput** d_reim, float** d_sum_re, float** d_sum_im,
	int** d_DivRayTubeNum, int** d_sum_gmem, int** d_sum_Gmem, int** d_squares_pred);

void FreeOnGPU(Direction* d_rays1, Square* d_squares1, Direction* d_rays2, Square* d_squares2, RayBeamInfo* d_effrays, Vector* d_center,
	Vector* d_axis, MatStruct* d_transMat, ReimOutput* d_reim, float* d_sum_re, float* d_sum_im, float* d_sum_sre, float* d_sum_sim, int* d_DivRayTubeNum, int* d_sum_gmem, int* d_sum_Gmem, int* d_squares_pred);

void free_data(Prim_Box *d_array, KD_Node_V *d_root, Element *d_points, Triangle *d_triangles);

void MemsetOnGPU(int dwidth, int dheight, Direction** d_rays1, Square** d_squares1, Direction** d_rays2, Square** d_squares2,
	RayBeamInfo** d_effrays, Vector** d_center, Vector** d_axis, MatStruct** d_transMat, ReimOutput** d_reim,
	float** d_sum_re, float** d_sum_im, float** d_sum_sre,float** d_sum_sim, int** d_DivRayTubeNum, int** d_sum_gmem, int** d_sum_Gmem, int** d_squares_pred);
