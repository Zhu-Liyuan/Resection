// Photogrammetry_Resection.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <Eigen/Dense>
#include <math.h>
using namespace Eigen;     //调用了矩阵运算库Eigen 
using namespace std;
double* Collinear_equation(double, double, double, double, double, double, double, double, double, double);
FILE* fp;
int main()
{
	//根据题目所给条件和数值创建变量和初始化
	double x[4] = { -0.08615,-0.05340,-0.01478,0.01046 };
	double y[4] = { -0.06899,0.08221,-0.07663,0.06443 };
	double cX[4] = { 36589.41,37631.08,39100.97,40426.54 }; //已知4个控制点X,Y,Z坐标
	double cY[4] = { 25273.32,31324.51,24934.98,30319.81 };
	double cZ[4] = { 2195.17,728.69,2386.50,757.31 };
	double x0 = 0, y0 = 0, f = 0.15324, m = 15000;
	double Xs = 0, Ys = 0, Zs = 0, q = 0, w = 0, k = 0, H = f * m;

	//计算投影中心的初值
	for (int i = 0; i < 4; i++) {
		Xs += cX[i];
		Ys += cY[i];
		Zs += cZ[i];
	}
	Xs /= 4; Ys /= 4; Zs = Zs / 4 + H;
	cout << "初始摄影中心坐标如下：" << endl;
	cout << Xs << "," << Ys << "," << Zs << endl;
	

	//A、V、X、L 矩阵初始化
	MatrixXd V(8, 1), A(8, 6), X(6, 1), L(8, 1), R(3, 3);

	int times;//传递次数
	

	for (int i = 0;; i++) {

		R << cos(q) * cos(k) - sin(q) * sin(w) * sin(k), -cos(q) * sin(k) - sin(q) * sin(w) * cos(k), -sin(q) * cos(w),
			cos(w)* sin(k), cos(w)* cos(k), -sin(w),
			sin(q)* cos(k) + cos(q) * sin(w) * sin(k), -sin(q) * sin(k) + cos(q) * sin(w) * cos(k), cos(q)* cos(w);
		
		for (int j = 0; j < 4; j++) {
			double X_ = R(0, 0) * (cX[j] - Xs) + R(1, 0) * (cY[j] - Ys) + R(2, 0) * (cZ[j] - Zs);
			double Y_ = R(0, 1) * (cX[j] - Xs) + R(1, 1) * (cY[j] - Ys) + R(2, 1) * (cZ[j] - Zs);
			double Z_ = R(0, 2) * (cX[j] - Xs) + R(1, 2) * (cY[j] - Ys) + R(2, 2) * (cZ[j] - Zs);
			double xx = x[j] - x0;
			double yy = y[j] - y0;
			A(2 * j, 0) = (R(0, 0) * f + R(0, 2) * xx) / Z_;
			A(2 * j, 1) = (R(1, 0) * f + R(1, 2) * xx) / Z_;
			A(2 * j, 2) = (R(2, 0) * f + R(2, 2) * xx) / Z_;
			A(2 * j, 3) = yy * sin(w) - (xx * (xx * cos(k) - yy * sin(k)) / f + f * cos(k)) * cos(w);
			A(2 * j, 4) = -f * sin(k) - xx * (xx * sin(k) + yy * cos(k)) / f;
			A(2 * j, 5) = yy;
			A(2 * j + 1, 0) = (R(0, 1) * f + R(0, 2) * yy) / Z_;
			A(2 * j + 1, 1) = (R(1, 1) * f + R(1, 2) * yy) / Z_;
			A(2 * j + 1, 2) = (R(2, 1) * f + R(2, 2) * yy) / Z_;
			A(2 * j + 1, 3) = -xx * sin(w) - (yy * (xx * cos(k) - yy * sin(k)) / f - f * sin(k)) * cos(w);
			A(2 * j + 1, 4) = -f * cos(k) - yy * (xx * sin(k) + yy * cos(k)) / f;
			A(2 * j + 1, 5) = -xx;
		}

		MatrixXd AT = A.transpose();
		
		double* temp;
		for (int i = 0; i < 4; i++) {
			temp = Collinear_equation(Xs, Ys, Zs, cX[i], cY[i], cZ[i], q, w, k, f);
			
			L(2 * i, 0) = x[i]-temp[0];
			L(2 * i + 1, 0) = y[i]-temp[1];
		}
		X = (AT * A).inverse() * AT * L;

		//更新迭代参数
		Xs += X(0, 0); Ys += X(1, 0); Zs += X(2, 0);
		q += X(3, 0); w += X(4, 0); k += X(5, 0);

		if ((X(3, 0)) < 0.000001 && (X(4, 0)) < 0.000001 && (X(5, 0)) < 0.000001 && i > 100) 
		{
			times = i;
			break;
		}

	}
	V = A * X - L;
	MatrixXd AT1 = A.transpose();
	MatrixXd Q = (AT1 * A).inverse();
	double m0 = sqrt((V.transpose() * V)(0, 0) / (2 * 4 - 6));
	cout << "\n解算之后的摄影中心坐标如下：" << endl;
	printf("Xs=%.2f\nYs=%.2f\nZs=%.2f\nq=%f\nw=%f\nk=%f\n", Xs, Ys, Zs, q, w, k);
	cout << "对应的旋转矩阵为：" << endl << R << endl;
	
	

	printf("\n\n单位权中误差为：%f\n\n", m0);
	printf("外方位元素的各精度为：\n%f\n%f\n%f\n%f\n%f\n%f\n", m0*sqrt(Q(0, 0)), m0*sqrt(Q(1, 1)), m0*sqrt(Q(2, 2)), m0*sqrt(Q(3, 3)), m0*sqrt(Q(4, 4)), m0*sqrt(Q(5, 5)));
	
	errno_t err;
	err = fopen_s(&fp, "后方交会结果输出.txt", "w+");
	fprintf_s(fp, "空间后方交会结果输出\n\n外方位元素解为:\nXs=%.2f\nYs=%.2f\nZs=%.2f\n", Xs, Ys, Zs);
	fprintf_s(fp, "\n外方位角元素构成的旋转矩阵的解R：\n%9f  %9f  %9f\n%9f  %9f  %9f\n%9f  %9f  %9f\n\n", R(0, 0), R(0, 1), R(0, 2), R(1, 0), R(1, 1), R(1, 2), R(2, 0), R(2, 1), R(2, 2));
	fprintf_s(fp, "迭代次数为:%d\n", times);
	fprintf_s(fp, "\n解算单位权中误差为:%f", m0);
	
	system("pause");
	return 0;
}


double* Collinear_equation(double Xs, double Ys, double Zs, double XA, double YA, double ZA, double  q, double w, double k, double f) 
{
	MatrixXd r(3, 3);
	r << cos(q) * cos(k) - sin(q) * sin(w) * sin(k), -cos(q) * sin(k) - sin(q) * sin(w) * cos(k), -sin(q) * cos(w),
		cos(w)* sin(k), cos(w)* cos(k), -sin(w),
		sin(q)* cos(k) + cos(q) * sin(w) * sin(k), -sin(q) * sin(k) + cos(q) * sin(w) * cos(k), cos(q)* cos(w);
	
	
	double x = -f * (r(0, 0) * (XA - Xs) + r(1, 0) * (YA - Ys) + r(2, 0) * (ZA - Zs)) / (r(0, 2) * (XA - Xs) + r(1, 2) * (YA - Ys) + r(2, 2) * (ZA - Zs));
	double y = -f * (r(0, 1) * (XA - Xs) + r(1, 1) * (YA - Ys) + r(2, 1) * (ZA - Zs)) / (r(0, 2) * (XA - Xs) + r(1, 2) * (YA - Ys) + r(2, 2) * (ZA - Zs));
	double* s = new double[2];
	s[0] = x;
	s[1] = y;
	return s;
}