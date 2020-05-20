//王瑞 2020.5.20
#ifndef __ELLIPSOID_REPAIR_H__
#define __ELLIPSOID_REPAIR_H__

typedef struct
{
	double x;
	double y;
	double z;
}s_ellipsoid_sample_point_t;//传感器读出的采样点



typedef struct
{

	double augmented_matrix[7][6];//6行7列的增广矩阵 [A|b]
	double array_solve[6];

}s_ellipsoid_repair_matrix_t;//增广矩阵和解向量

typedef struct
{
	//offset
	double X;
	double Y;
	double Z;

	//scale coefficient
	double a;
	double b;
	double c;
}s_ellipsoid_param_t;//椭球校正计算完成后的参数


int ellipsoid_repair(s_ellipsoid_sample_point_t* sp_sample_points, int N, s_ellipsoid_param_t* sp_repair_param);


#endif

