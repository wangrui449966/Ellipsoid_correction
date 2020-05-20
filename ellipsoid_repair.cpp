// 王瑞  2020.5.20

#include "ellipsoid_repair.h"
#include <stdint.h>
#include <math.h>


#include <iostream>
#include <stdio.h>
using namespace std;
static void debug_print_matrix(double mat[7][6])
{
	int i, j;
	for (j = 0; j < 6; j++)
	{
		for (i = 0; i < 7; i++)
		{
			if (fabs(mat[i][j]) > 0.001)
			{
				printf("%.3lf\t", mat[i][j]);
			}
			else
			{
				printf("0\t");
			}
		}
		cout << endl;
	}
	cout << endl << endl;
}

/*
@param sp_sample_points:采样点数组首地址
@param N:采样点个数，必须有6个不同的点来保证系数矩阵满秩
@param sp_repair_param:计算结束返回的椭球校正参数
*/
int ellipsoid_repair(s_ellipsoid_sample_point_t* sp_sample_points, int N,s_ellipsoid_param_t *sp_repair_param)
{
	s_ellipsoid_repair_matrix_t s_matrix;
	int i, j;

	if (N < 6)
	{
		return -1;//N小于6
	}

	double pre_norm_coefficient;//用第一个采样点来计算的预归一化系数，让模值接近与1，防止计算过程中导致数据太大或者太小，最终的时候再除以这个系数
	pre_norm_coefficient = 1.0 / sqrt(pow(sp_sample_points[0].x, 2) + pow(sp_sample_points[0].y, 2) + pow(sp_sample_points[0].z, 2));
	
	//初始化s_matrix
	for (i = 0; i < sizeof(s_ellipsoid_repair_matrix_t); i++)
	{
		((uint8_t*)&s_matrix)[i] = 0;
	}

	//增广矩阵赋值
	for (i = 0; i < N; i++)
	{
		double xi, yi, zi;
		xi = sp_sample_points[i].x * pre_norm_coefficient;
		yi = sp_sample_points[i].y* pre_norm_coefficient;
		zi = sp_sample_points[i].z * pre_norm_coefficient;
		
		s_matrix.augmented_matrix[0][0] += pow(yi, 4);
		s_matrix.augmented_matrix[1][0] += pow(yi, 2) * pow(zi, 2);
		s_matrix.augmented_matrix[2][0] += pow(xi, 1) * pow(yi, 2);
		s_matrix.augmented_matrix[3][0] += pow(yi, 3);
		s_matrix.augmented_matrix[4][0] += pow(yi, 2) * pow(zi, 1);
		s_matrix.augmented_matrix[5][0] += pow(yi, 2);
		s_matrix.augmented_matrix[6][0] += -pow(xi, 2) * pow(yi, 2);


		s_matrix.augmented_matrix[0][1] += pow(yi, 2) * pow(zi, 2);
		s_matrix.augmented_matrix[1][1] += pow(zi, 4);
		s_matrix.augmented_matrix[2][1] += pow(xi, 1) * pow(zi, 2);
		s_matrix.augmented_matrix[3][1] += pow(yi, 1) * pow(zi, 2);
		s_matrix.augmented_matrix[4][1] += pow(zi, 3);
		s_matrix.augmented_matrix[5][1] += pow(zi, 2);
		s_matrix.augmented_matrix[6][1] += -pow(xi, 2) * pow(zi, 2);


		s_matrix.augmented_matrix[0][2] += pow(xi, 1) * pow(yi, 2);
		s_matrix.augmented_matrix[1][2] += pow(xi, 1) * pow(zi, 2);
		s_matrix.augmented_matrix[2][2] += pow(xi, 2);
		s_matrix.augmented_matrix[3][2] += pow(xi, 1) * pow(yi, 1);
		s_matrix.augmented_matrix[4][2] += pow(xi, 1) * pow(zi, 1);
		s_matrix.augmented_matrix[5][2] += pow(xi, 1);
		s_matrix.augmented_matrix[6][2] += -pow(xi, 3);

		s_matrix.augmented_matrix[0][3] += pow(yi, 3);
		s_matrix.augmented_matrix[1][3] += pow(yi, 1) * pow(zi, 2);
		s_matrix.augmented_matrix[2][3] += pow(xi, 1) * pow(yi, 1);
		s_matrix.augmented_matrix[3][3] += pow(yi, 2);
		s_matrix.augmented_matrix[4][3] += pow(yi, 1) * pow(zi, 1);
		s_matrix.augmented_matrix[5][3] += pow(yi, 1);
		s_matrix.augmented_matrix[6][3] += -pow(xi, 2) * pow(yi, 1);

		s_matrix.augmented_matrix[0][4] += pow(yi, 2) * pow(zi, 1);
		s_matrix.augmented_matrix[1][4] += pow(zi, 3);
		s_matrix.augmented_matrix[2][4] += pow(xi, 1) * pow(zi, 1);
		s_matrix.augmented_matrix[3][4] += pow(yi, 1) * pow(zi, 1);
		s_matrix.augmented_matrix[4][4] += pow(zi, 2);
		s_matrix.augmented_matrix[5][4] += pow(zi, 1);
		s_matrix.augmented_matrix[6][4] += -pow(xi, 2) * pow(zi, 1);


		s_matrix.augmented_matrix[0][5] += pow(yi, 2);
		s_matrix.augmented_matrix[1][5] += pow(zi, 2);
		s_matrix.augmented_matrix[2][5] += pow(xi, 1);
		s_matrix.augmented_matrix[3][5] += pow(yi, 1);
		s_matrix.augmented_matrix[4][5] += pow(zi, 1);
		s_matrix.augmented_matrix[5][5] += 1;
		s_matrix.augmented_matrix[6][5] += -pow(xi, 2);
	}

	cout << "增广矩阵" << endl;
	debug_print_matrix(s_matrix.augmented_matrix);

	//化为行最简阶梯形矩阵
	//先化上三角
	for (i = 0; i < 5; i++)
	{
		for (j = i + 1; j < 6; j++)
		{
			double k;//某一行乘k倍加入到另一行
			k = -s_matrix.augmented_matrix[i][j] / s_matrix.augmented_matrix[i][i];
			int L;
			for (L = 0; L < 7; L++)
			{
				s_matrix.augmented_matrix[L][j] += k * s_matrix.augmented_matrix[L][i];
			}
		}
	}
	cout << "上三角" << endl;
	debug_print_matrix(s_matrix.augmented_matrix);
	//再化为行最简阶梯
	for (i = 5; i > 0; i--)
	{
		for (j = i - 1; j >= 0; j--)
		{
			double k;////某一行乘k倍加入到另一行
			k = -s_matrix.augmented_matrix[i][j] / s_matrix.augmented_matrix[i][i];
			int L;
			for (L = 0; L < 7; L++)
			{
				s_matrix.augmented_matrix[L][j] += k * s_matrix.augmented_matrix[L][i];
			}
		}
	}
	cout << "行最简" << endl;
	debug_print_matrix(s_matrix.augmented_matrix);

	//解向量
	for (i = 0; i < 6; i++)
	{
		s_matrix.array_solve[i] = s_matrix.augmented_matrix[6][i] / s_matrix.augmented_matrix[i][i];
	}

	double A, B, C, D, E, F;
	double X, Y, Z, a, b, c;
	A = s_matrix.array_solve[0];
	B = s_matrix.array_solve[1];
	C = s_matrix.array_solve[2];
	D = s_matrix.array_solve[3];
	E = s_matrix.array_solve[4];
	F = s_matrix.array_solve[5];

	if (A <= 0 || B <= 0)
	{
		return -2;//A和B都应该大于0
	}



	X = -0.5 * C;
	Y = D / (-2.0 * A);
	Z = -E / (2.0 * B);
	a = sqrt(X * X + A * Y * Y + B * Z * Z - F);
	b = sqrt(a * a / A);
	c = sqrt(a * a / B);

	X /= pre_norm_coefficient;
	Y /= pre_norm_coefficient;
	Z /= pre_norm_coefficient;
	a /= pre_norm_coefficient;
	b /= pre_norm_coefficient;
	c /= pre_norm_coefficient;

	sp_repair_param->a = a;
	sp_repair_param->b = b;
	sp_repair_param->c = c;
	sp_repair_param->X = X;
	sp_repair_param->Y = Y;
	sp_repair_param->Z = Z;


	cout << "X= " << X << endl;
	cout << "Y= " << Y << endl;
	cout << "Z= " << Z << endl;
	cout << "a= " << a << endl;
	cout << "b= " << b << endl;
	cout << "c= " << c << endl;

	return 0;
}
