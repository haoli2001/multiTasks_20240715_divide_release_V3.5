#ifndef REFLECTCOEFF_H
#define REFLECTCOEFF_H
#define A2R 3.1415926535897/180
/**************************
名称：ReflectCoeff_2.h
描述：反射系数函数头文件
***************************/
struct comp1
{
  float re;
  float im;
};

__host__ float ReflectCoeff_2(float f, float theta);

#endif