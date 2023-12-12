#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <iomanip>

using namespace std;

double func1(double x1, double x2)
{
	return log(1 + (x1 + x2) / 5.0) - sin(x2 / 3.0) + 1.1;
}

double func2(double x1, double x2)
{
	return cos(x1 * x2 / 6.0) - x2 + 0.5;
}

double der00(double x1, double x2)
{
	return (1 / (x1 + x2 + 5.0));
}

double der01(double x1, double x2)
{
	return (1.0 / (x1 + x2 + 5.0) - (1.0 / 3.0) * cos(x2 / 3.0));
}

double der10(double x1, double x2)
{
	return (-1 * (1.0 / 6.0) * x2 * sin(x1 * x2 / 6.0));
}

double der11(double x1, double x2)
{
	return (-1 * (1.0 / 6.0) * x1 * sin(x1 * x2 / 6.0) - 1);
}

double* methodG(double** full, int n, int m)
{
	//деление построчно
	double elem, max;
	int coordStr;

	for (int i = 0; i < n; i++)
	{
		max = 0;
		coordStr = 0;

		for (int j = i; j < n; j++) //нахождение максимума
		{
			if (abs(full[j][i]) > max)
			{
				max = abs(full[j][i]);
				coordStr = j;
			}
		}

		if (max == 0)
		{
			cout << "Max in Gauss method == 0";
			exit;
		}

		if (max > abs(full[i][i])) //меняем местами 1 строку и строку с наибольшим 1 элементом
		{
			double* ptr = full[i];
			full[i] = full[coordStr];
			full[coordStr] = ptr;
		}

		elem = full[i][i];

		for (int j = i; j < m; j++)
		{
			full[i][j] /= elem;   //делим строку на элемент
		}

		for (int j = i + 1; j < n; j++) //делим следующие строки
		{
			elem = full[j][i];

			for (int k = i; k < m; k++)
			{
				full[j][k] -= elem * full[i][k];
			}
		}
	}

	//поиск решения

	double* ans = new double[n];
	ans[n - 1] = full[n - 1][m - 1]; //

	for (int i = n - 2; i >= 0; i--)
	{
		ans[i] = full[i][n];

		for (int j = i + 1; j < m - 1; j++)
		{
			ans[i] -= full[i][j] * ans[j];
		}
	}

	return ans;
}

// double** jacobFunc(double x1, double x2)
// {
// 	double Jacob_func1[2][2]{ {der00(x1 ,x2), der01(x1 ,x2)}, {der10(x1 ,x2), der11(x1 ,x2)} };
// 	double** Jacob_func = new double* [2];

// 	for (int i = 0; i < 2; i++)
// 	{
// 		Jacob_func[i] = new double[2];
// 		for (int j = 0; j < 2; j++)
// 		{
// 			Jacob_func[i][j] = Jacob_func1[i][j];
// 		}
// 	}
	
// 	return Jacob_func;
// }

double** jacobFunc(double x1, double x2) 
{
	double Jacob_func1[2][2] {{0}};
	double** Jacob_func = new double* [2];
	double M = 0.01;

	Jacob_func1[0][0] = (func1(x1 + M, x2) - func1(x1, x2)) / M;
	Jacob_func1[0][1] = (func1(x1, x2 + M) - func1(x1, x2)) / M;
	Jacob_func1[1][0] = (func2(x1 + M, x2) - func2(x1, x2)) / M;
	Jacob_func1[1][1] = (func2(x1, x2 + M) - func2(x1, x2)) / M;

	for (int i = 0; i < 2; i++)
	{
		Jacob_func[i] = new double[2];
		for (int j = 0; j < 2; j++)
		{
			Jacob_func[i][j] = Jacob_func1[i][j];
		}
	}

	return Jacob_func;
}

double** forG(double** Jacob_func, double* results)
{
	double** G = new double* [2];

	for (int i = 0; i < 2; i++)
	{
		G[i] = new double[3];
		G[i][2] = -1 * results[i];
		for (int j = 0; j < 2; j++)
		{
			G[i][j] = Jacob_func[i][j];
		}
	}

	return G;
}

void out(int iter, double x1, double x2, double temp1, double temp2)
{
	double max1 = 0, max2 = 0;
	if (func1(temp1, temp2) > max1)
	{
		max1 = func1(temp1, temp2);
	}
	else if (func2(temp1, temp2) > max1)
	{
		max1 = func2(temp1, temp2);
	}

	if (abs(x1) < 1.0)
	{
		if ((x1 - temp1) > max2)
		{
			max2 = x1 - temp1;
		}
		else if ((x2 - temp2) > max2)
		{
			max2 = x2 - temp2;
		}
	}
	else
	{
		if ((x1 - temp1)/x1 > max2)
		{
			max2 = (x1 - temp1)/x1;
		}
		else if ((x2 - temp2)/x2 > max2)
		{
			max2 = (x2 - temp2)/x2;
		}
	}

	cout << "\t" << iter << "\t" << setprecision(9)<< max1 << "\t"<< setprecision(9) << max2 << "\n\n";
}

int Max_func(int x, int y){
	double max_func;
	double f1 = func1(x, y);
	double f2 = func2(x, y);

	abs(f1 > f2) ?
		max_func = f1 :
		max_func = f2;

	return max_func;
}

int Max_ans(int x, int y){
	double max_ans;
	abs(x > y) ?
			max_ans = x :
			max_ans = y;
	
	return max_ans;
}

int main()
{

	double eplission = pow(10, -9);
	double x1 = 1.0, x2 = 1.0, temp1 = 0, temp2 = 0;
	int iter = 1;
	cout << "\tI\td1\td2\n\n";
	while (true)
	{
		temp1 = x1;
		temp2 = x2;
		double** Jacob_func = jacobFunc(x1, x2);
		double equation_results[2]{ func1(x1 ,x2), func2(x1 ,x2) };
		double* ans = methodG(forG(Jacob_func, equation_results), 2, 3);
		x1 += ans[0];
		x2 += ans[1];
		out(iter, x1, x2, temp1, temp2);
		
		if (Max_ans(ans[0], ans[1]) < eplission && Max_func(ans[0], ans[1]) < eplission)
		{
			break;
		}
		if (iter > 50)
		{
			cout << "Iterations error";
			break;
		}
		iter++;
	}
	cout << "\tx1 = \t" << x1 << "\tx2 = \t" << x2;
}
