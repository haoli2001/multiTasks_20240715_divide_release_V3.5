#ifndef REFLECTCOEFF_H
#define REFLECTCOEFF_H
#define A2R 3.1415926535897/180
/**************************
���ƣ�ReflectCoeff_2.h
����������ϵ������ͷ�ļ�
***************************/
struct comp1
{
  float re;
  float im;
};

__host__ float ReflectCoeff_2(float f, float theta);

#endif