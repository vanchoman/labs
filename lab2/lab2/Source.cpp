#include <iostream>
#include <vector>

using namespace std;

bool CheckProportionality(vector<vector<double>> A, vector<double> b);
bool CheckRowProportionality(vector<double> firstRows, vector<double> secondRows);
vector<double> SolveGaussMethod(vector<vector<double>> A, vector<double> b);
vector<vector<double>> AnalyticalMethodForJacobiMatrix(double x1, double x2);
vector<vector<double>> CourseDifferenceMethodForJacobiMatrix(double x1, double x2, double relIncrement);
void SolveNewtonMethod(double x1, double x2, const double firstSolutionError, const double secondSolutionError,
	const int maxNumberIterations, double relativeIncrement);
double f1(double x1, double x2) { return sin(x1 + 1) - x2 - 1; }
double f1dx1(double x1, double x2) { return cos(x1 + 1); }
double f1dx2(double x1, double x2) { return -1.0; }
double f2(double x1, double x2) { return 2 * x1 + cos(x2) - 2; }
double f2dx1(double x1, double x2) { return 2; }
double f2dx2(double x1, double x2) { return -sin(x2); }


bool CheckProportionality(vector<vector<double>> A, vector<double> b)
{
	bool isProportional = false;
	int vectorSize = A.size();
	vector<vector<double>> A_b(vectorSize, vector<double>(vectorSize + 1, 0));
	for (int i = 0; i < vectorSize; i++)
	{
		for (int j = 0; j < vectorSize + 1; j++)
		{
			if (j == vectorSize)
				A_b[i][j] = b[i];
			else
				A_b[i][j] = A[i][j];
		}
	}
	for (int i = 0; i < vectorSize - 1; i++)
	{
		for (int j = i + 1; j < vectorSize; j++)
		{
			isProportional = CheckRowProportionality(A_b[i], A_b[j]);
			if (isProportional)
				return true;
		}
	}
	return false;
}

bool CheckRowProportionality(vector<double> firstRows, vector<double> secondRows)
{
	int rowSize = firstRows.size();
	double proportionalityFactor = firstRows[0] / secondRows[0];
	for (int i = 1; i < rowSize; i++)
	{
		if (proportionalityFactor == firstRows[i] / secondRows[i])
			continue;
		else
			return false;
	}
	return true;
}

vector<double> SolveGaussMethod(vector<vector<double>> A, vector<double> b)
{
	if (CheckProportionality(A, b))
	{
		cout << "The rows in the matrix are proportional!";
		exit(1);
	}
	int vectorSize = A.size();
	for (int i = 0; i < vectorSize; i++)
	{
		int k = i;
		for (int j = i + 1; j < vectorSize; j++)
		{
			if (abs(A[j][i]) > abs(A[k][i]))
				k = j;
		}
		swap(A[i], A[k]);
		swap(b[i], b[k]);
		double div = A[i][i];
		for (int j = i; j < vectorSize; j++)
			A[i][j] /= div;
		b[i] /= div;
		for (int j = i + 1; j < vectorSize; j++)
		{
			double mult = A[j][i];
			for (int k = i; k < vectorSize; k++)
				A[j][k] -= mult * A[i][k];
			b[j] -= mult * b[i];
		}
	}
	vector<double> vectorOfRoots(vectorSize);
	for (int i = vectorSize - 1; i >= 0; i--)
	{
		vectorOfRoots[i] = b[i];
		for (int j = i + 1; j < vectorSize; j++)
			vectorOfRoots[i] -= A[i][j] * vectorOfRoots[j];
	}
	return vectorOfRoots;
}

vector<vector<double>> AnalyticalMethodForJacobiMatrix(double x1, double x2)
{
	return { {f1dx1(x1, x2), f1dx2(x1, x2)}, { f2dx1(x1, x2), f2dx2(x1, x2) } };
}

vector<vector<double>> CourseDifferenceMethodForJacobiMatrix(double x1, double x2, double relIncrement)
{
	return vector<vector<double>> { {(f1(x1 + x1 * relIncrement, x2) - f1(x1, x2)) / relIncrement / x1, (f1(x1, x2 + x2 * relIncrement) - f1(x1, x2)) / relIncrement / x2},
		{ (f2(x1 + x1 * relIncrement, x2) - f2(x1, x2)) / relIncrement / x1, (f2(x1, x2 + x2 * relIncrement) - f2(x1, x2)) / relIncrement / x2 } };
}

void SolveNewtonMethod(double x1, double x2, const double firstSolutionError, const double secondSolutionError,
	const int maxNumberIterations, double relativeIncrement)
{
	double firstDelta = max(abs(f1(x1, x2)), abs(f2(x1, x2)));
	double secondDelta = 1;
	if (relativeIncrement != NULL)
		cout << "Relative increment: " << relativeIncrement << ";\n\n";
	int iteration = 0;
	while ((firstDelta > firstSolutionError || secondDelta > secondSolutionError) && iteration < maxNumberIterations)
	{
		iteration++;
		printf("%d: delta of x1: %.*f; delta of x2: %.*f;\n", iteration, 10, firstDelta, 10, secondDelta);
		vector<double> residualVector{ -f1(x1, x2), -f2(x1, x2) };
		vector<vector<double>> JacobiMatrix;
		if (relativeIncrement == NULL)
			JacobiMatrix = AnalyticalMethodForJacobiMatrix(x1, x2);
		else
			JacobiMatrix = CourseDifferenceMethodForJacobiMatrix(x1, x2, relativeIncrement);
		vector<double> vectorOfSolution = SolveGaussMethod(JacobiMatrix, residualVector);
		x1 += vectorOfSolution[0];
		x2 += vectorOfSolution[1];
		firstDelta = abs(residualVector[0]);
		for (int i = 1; i < residualVector.size(); i++)
		{
			if (firstDelta < abs(residualVector[i]))
				firstDelta = abs(residualVector[i]);
		}
		double firstMax = abs(x1) < 1 ? abs(vectorOfSolution[0]) : abs(vectorOfSolution[0] / x1);
		double secondMax = abs(x2) < 1 ? abs(vectorOfSolution[1]) : abs(vectorOfSolution[1] / x2);
		secondDelta = max(firstMax, secondMax);
	}
	printf("\nnumber of iteration: %d. \nfirst x on this iteration: %.*f; \nsecond x on this iteration: %.*f.\n",
		iteration, 15, x1, 15, x2);
	cout << "\n=====================================================\n\n";
}

int main()
{
	const int maxNumberIterations = 100;
	const double firstSolutionError = pow(10, -9), secondSolutionError = pow(10, -9);
	double x1 = 1, x2 = 1;
	cout << "General information for all solutions: \n";
	printf("Starting point of approximation for all solutions of the system: (%.*f; %.*f);\n", 0, x1, 0, x2);
	printf("Solution error for x1: %.*f;\nSolution error for x2: %.*f;\n",
		9, firstSolutionError, 9, secondSolutionError);
	cout << "Max number of iterations: " << maxNumberIterations <<
		".\n\n=====================================================\n\n";
	vector<double> vectorOfRelativeIncrement{ 0.01, 0.05, 0.1 };
	SolveNewtonMethod(x1, x2, firstSolutionError, secondSolutionError, maxNumberIterations, NULL);
	for (int i = 0; i < vectorOfRelativeIncrement.size(); i++)
		SolveNewtonMethod(x1, x2, firstSolutionError, secondSolutionError, maxNumberIterations,
			vectorOfRelativeIncrement[i]);
	return 0;
}