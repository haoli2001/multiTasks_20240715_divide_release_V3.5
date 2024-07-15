#include "ReflectCoeff_2.h"

__device__ __host__ comp1 cdiv(comp1 z1,comp1 z2)       
{ double x1,x2,y1,y2;
  comp1 z;
  x1=z1.re;
  x2=z2.re;
  y1=z1.im;
  y2=z2.im;
  z.re=(x1*x2+y1*y2)/(x2*x2+y2*y2);
  z.im=(x2*y1-y2*x1)/(x2*x2+y2*y2);
  return z;
}

__device__ __host__ comp1 cmul(comp1 z1,comp1 z2)
{ double x1,x2,y1,y2;
  comp1 z;
  x1=z1.re;
  x2=z2.re;
  y1=z1.im;
  y2=z2.im;
  z.re=x1*x2-y1*y2;
  z.im=x1*y2+y1*x2;
  return z;
}

__device__ __host__ double cabs(comp1 z)                      
{
	double x,y;
	x=z.re;
	y=z.im;
	return sqrt(x*x+y*y);
}

//����˫������T��ķ���ϵ��
__host__ float ReflectCoeff_2(float f, float theta)
{
	double theta0, theta1, theta2;
	double c0, c1, c2, c3;
	double rou0, rou1, rou2, rou3;
	double Ee2;
	double Z0, Z3, Z4;
	double eta1, eta2;
	double k1, k2, k3;
	double d1, d2, d3;
	double h1, h2, h3;
	double R, AR;
	double phi;
	double AR0, AR1;
	comp1 temp1, temp2, Z23, Z12, Z01, Zin, cw1, cw2, Z1, Z2;

	theta0 = theta;  //����Ƕ�
	rou0 = 1000;     //ˮ�ܶȣ���λkg/m^3��
	c0 = 1500;       //ˮ�����٣���λm/s��
	Z0 = rou0 * c0;  //ˮ���迹

	rou1 = 1039;     //�𽺣�����1���ܶ�
	c1 = 1470;     //����1�еĵ�Ч����
	eta1 = 0.4;      //�������
	phi = atan(eta1);
	cw1.re = c1 * pow((1 * 1 + eta1*eta1), 0.25)*cos(phi / 2);
	cw1.im = c1 * pow((1 * 1 + eta1*eta1), 0.25)*sin(phi / 2);//ճ�����ϲ���
	Z1.re = rou1 * cw1.re;    //���ϵ��ݲ����迹
	Z1.im = rou1 * cw1.im;    //���ϵ��ݲ����迹
	k1 = 2 * 180 * A2R * f / c1;//����1����

	rou2 = 1090;     //�𽺣�����1���ܶ�
	eta2 = 0.5;      //�������
	phi = atan(eta2);
	Ee2 = 1e9;      //����2������ģ��
	c2 = sqrt(Ee2 / rou2);//����2�еĵ�Ч����
	cw2.re = c2 * pow((1 * 1 + eta2*eta2), 0.25)*cos(phi / 2);
	cw2.im = c2 * pow((1 * 1 + eta2*eta2), 0.25)*sin(phi / 2);//ճ�����ϲ���
	Z2.re = rou2 * cw2.re;    //���ϵ��ݲ����迹
	Z2.im = rou2 * cw2.im;    //���ϵ��ݲ����迹
	k2 = 2 * 180 * A2R * f / c2;//���ʲ���

	rou3 = 7850;     //�֣��ĵף��ܶ�
	c3 = 5200;       //����2�еĵ�Ч����
	Z3 = rou3 * c3;  //�ĵ�(�ֵ��迹����ˮ��20��)
	k3 = 2 * 180 * A2R * f / c3;//���ʲ���,�ĵ�Ϊ�ְ�ʱ�Ĵ����ٶ�

	d1 = 2e-3;       //����1���
	d2 = 2e-3;       //����2���
	d3 = 3e-3;       //�ĵ׺��

					 //////����Snell���䶨�����������ʲ��е�����Ƕ�
	AR0 = sin(theta0) * c0 / c1;  //����ȫ����ϵ��
	if (fabs(AR0) > 1)
	{
		theta1 = 90 * A2R;
	}
	else
	{
		theta1 = asin(sin(theta0) * c0 / c1);
	}

	AR1 = sin(theta1) * c1 / c2;  //����ȫ����ϵ��
	if (fabs(AR1) > 1)
	{
		theta2 = 90 * A2R;
	}
	else
	{
		theta2 = asin(sin(theta1) * c1 / c2);
	}
	h1 = tan(k1 * cos(theta0) * d1);
	h2 = tan(k2 * cos(theta1) * d2);
	h3 = tan(k3 * cos(theta2) * d3);
	Z4 = Z0;

	temp1.im = Z3 * h3;
	temp1.re = Z4;
	temp2.im = Z4 * h3;
	temp2.re = Z3;
	Z23 = cdiv(temp1, temp2);
	Z23.re = Z3 * Z23.re;
	Z23.im = Z3 * Z23.im;


	temp1.im = Z23.im + Z2.re * h2;
	temp1.re = Z23.re - h2*Z2.im;
	temp2.im = Z23.re * h2 + Z2.im;
	temp2.re = Z2.re - Z23.im * h2;
	Z12 = cdiv(temp1, temp2);
	Z12 = cmul(Z12, Z2);

	temp1.im = Z12.im + Z1.re * h1;
	temp1.re = Z12.re - h1*Z1.im;
	temp2.im = Z1.im + Z12.re * h1;
	temp2.re = Z1.re - Z12.im * h1;
	Z01 = cdiv(temp1, temp2);
	Z01 = cmul(Z01, Z1);
	Zin = Z01;
	//Zin = Z1 * (Z3 * (Z2 - h1 * h2 * Z1) + jay * Z2 * (Z2 * h2 + h1 * Z1))/(Z2 * (Z1 - h1 * h2 * Z2) + jay * Z3 * (Z1 * h2 + h1 * Z2));  //�����迹
	R = fabs((cabs(Zin) - Z0) / (cabs(Zin) + Z0)); //����ϵ��
	return R;
}
