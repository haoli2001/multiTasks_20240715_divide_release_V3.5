#ifndef HANDLERROR_H

#define HANDLERROR_H



#include "cuda_runtime.h"
#include <iostream>

#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

using namespace std;

static void HandleError(cudaError_t err, const char *file, int line)

{

	if (err != cudaSuccess)

	{

		cout << "在" << file << "的第" << line << "行出现错误！" << endl;

		cout << "错误代码是：" << cudaGetErrorString(err) << endl;

	}

}





#endif
