#include "pch.h"
#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;

double fun(double x);

double funk(int i, int j, int n);

double trapec(double a, double b, double(*f)(double), double e);

double simpson2(double a, double b, double(*f)(double), double e);

double simpson3(double a, double A, double b, double B, double(*f)(int, int, int), double e);

int main()
{
	setlocale(LC_ALL, "Russian");
	double a = 1.2, b = 2.471, a1 = -1, b1 = -1, A1 = 1, B1 = 1;
	double e = pow(10, -4);
	cout << "e = 10^-4" << endl;
	cout << "Трапеции: "; trapec(a, b, fun, e);
	cout << "	Симпсоны2: "; simpson2(a, b, fun, e);
	cout << "	Симпсоны3: "; simpson3(a1, A1, b1, B1, funk, e);
	e = pow(10, -5);
	cout << endl << "e = 10^-5" << endl;
	cout << "Трапеции: "; trapec(a, b, fun, e);
	cout << "	Симпсоны2: "; simpson2(a, b, fun, e);
	cout << "	Симпсоны3: "; simpson3(a1, A1, b1, B1, funk, e);
	cout << endl;
}

double fun(double x)
{
	double F = pow((1 + 2 * pow(x, 3)), 0.5);
	return F;
}

double funk(int i, int j, int n)
{
	double a = -1, b = -1, A = 1, B = 1;
	double h = (A - a) / (2 * n), k = (B - b) / (2 * n);
	double x = a + i * h, y = b + j * k;
	double F = 4 - pow(x, 2) - pow(y, 2);
	if (x > 1 || y > 1)
		return 0;
	else
		return F;
}

double trapec(double a, double b, double(*f)(double), double e)
{
	int n = 2;
	double h = (b - a) / n, I = f(a) + f(b), I2, x = a;
	for (int i = 0; i < n - 1; i++)
	{
		x += h;
		I += 2 * f(x);
	}
	I *= h / 2;

	do
	{
		n *= 2;
		h = (b - a) / n;
		I2 = I;
		I = f(a) + f(b);
		x = a;
		for (int i = 0; i < n - 1; i++)
		{
			x += h;
			I += 2 * f(x);
		}
		I *= h / 2;
	} while (abs(I - I2) > 3 * e);
	double R;

	R = (I2 - I) / (pow(0.5, 2) - 1);

	cout << I << " (погрешность " << R << ")";
	return I;
}

double simpson2(double a, double b, double(*f)(double), double e)
{
	int n = 4;
	double h = (b - a) / n, I = f(a) + f(b), I2, x = a;
	for (int i = 0; i < n - 1; i += 2)
	{
		x += h;
		I += 4 * f(x);
		x += h;
		I += 2 * f(x);
	}
	I *= h / 3;

	do
	{
		n *= 2;
		cout << n << endl;
		h = (b - a) / n;
		I2 = I;
		I = f(a) + f(b);
		x = a;
		for (int i = 0; i < n - 1; i += 2)
		{
			x += h;
			I += 4 * f(x);
			x += h;
			I += 2 * f(x);
		}
		I *= h / 3;
	} while (abs(I - I2) > 15 * e);
	double R;

	R = (I2 - I) / (pow(0.5, 4) - 1);

	cout << I << " (погрешность " << R << ")";
	return I;
}

double simpson3(double a, double A, double b, double B, double (*f)(int, int, int), double e)
{
	int n = 4, m = n;
	double h = (A - a) / (2 * n), k = (B - b) / (2 * m);
	double I = 0, I2;

	for (int i = 0; i <= n; i++)
	{
		for (int j = 0; j <= m; j++)
		{
			I += f(2 * i, 2 * j, n) + f(2 * i + 2, 2 * j, n) + f(2 * i + 2, 2 * j + 2, n) + f(2 * i, 2 * j + 2, n) +
				4 * (f(2 * i + 1, 2 * j, n) + f(2 * i + 2, 2 * j + 1, n) + f(2 * i + 1, 2 * j + 2, n) + f(2 * i, 2 * j + 1, n)) +
				16 * f(2 * i + 1, 2 * j + 1, n);
		}
	}
	I *= h * k / 9;

	do
	{
		n *= 2;
		h = (A - a) / (2 * n);
		k = (B - b) / (2 * m);
		I2 = I;
		I = 0;
		for (int i = 0; i <= n; i++)
		{
			for (int j = 0; j <= m; j++)
			{
				I += f(2 * i, 2 * j, n) + f(2 * i + 2, 2 * j, n) + f(2 * i + 2, 2 * j + 2, n) + f(2 * i, 2 * j + 2, n) +
					4 * (f(2 * i + 1, 2 * j, n) + f(2 * i + 2, 2 * j + 1, n) + f(2 * i + 1, 2 * j + 2, n) + f(2 * i, 2 * j + 1, n)) +
					16 * f(2 * i + 1, 2 * j + 1, n);
			}
		}
		I *= h * k / 9;
	} while (abs(I - I2) > 15 * e);

	cout << I;
	return I;
}
