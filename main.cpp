#include <iostream>
#include <stdlib.h>
#include <Windows.h>
#include <stdio.h>
#include <conio.h>
#include <stdint.h>
#include "ellipsoid_repair.h"
using namespace std;


void sample_points_set(s_ellipsoid_sample_point_t sample_poinsts[], int i, double x, double y, double z)
{
	sample_poinsts[i].x = x;
	sample_poinsts[i].y = y;
	sample_poinsts[i].z = z;
	printf("加入采样点：\t%.1lf\t%.1lf\t%.1lf\n", x, y, z);
}
void main()
{
	s_ellipsoid_sample_point_t sample_poinsts[6];
	s_ellipsoid_param_t s_repair_param;
	int N=0;
	if (1)
	{
		sample_points_set(sample_poinsts, N++, 1.3, -0.1, 0.2);
		sample_points_set(sample_poinsts, N++, -1.1, -0.1, 0.2);
		sample_points_set(sample_poinsts, N++, 0.1, 1.3, 0.2);
		sample_points_set(sample_poinsts, N++, 0.1, -1.5, 0.2);
		sample_points_set(sample_poinsts, N++, 0.1, -0.1, 1.1);
		sample_points_set(sample_poinsts, N++, 0.1, -0.1, -0.7);
	}
	cout << endl;
	ellipsoid_repair(sample_poinsts, 6, &s_repair_param);


}