#ifndef __SCALFUC_GPU_H__
#define __SCALFUC_GPU_H__
#include "common_struct.h"
#include "integral_gpu.h"

void scalfuc(RayBeamInfo* rays, int raysBeamNum, int ig, float* d_sum_sre, float* d_sum_sim, ConfigStruct config, Vector* d_center, Axis_slx New_receive_points);

#endif
