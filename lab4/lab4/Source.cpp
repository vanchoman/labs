#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>

using namespace std;

double MultiplyParametersOfVectors(vector<double> firstVector, vector<double> secondVector) // Σ(fV[i] * sV[i])
{
	double resultOfMult = 0;
	for (int i = 0; i < firstVector.size(); i++)
		resultOfMult += firstVector[i] * secondVector[i];
	return resultOfMult;
}

double SumVectorElements(vector <double> vectorOfElements)
{
	double sumOfVectorElements = 0;
	for (int i = 0; i < vectorOfElements.size(); i++)
		sumOfVectorElements += vectorOfElements[i];
	return sumOfVectorElements;
}

vector<double> RaiseElementsOfVectorToPow(vector<double> vectorOfElements, int degree)
{
	vector<double> vectorOfResults(vectorOfElements.size());
	for (int i = 0; i < vectorOfElements.size(); i++)
		vectorOfResults[i] = pow(vectorOfElements[i], degree);
	return vectorOfResults;
}

vector<double> SolveGaussMethod(vector<vector<double>> A, vector<double> b)
{
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

vector <double> SolveLSMApproximation(vector <double> experimentalValuesOfX, vector <double> experimentalValuesOfY)
{
	int numberOfExperimentalPoints = experimentalValuesOfX.size();
	int polynomialDegree = 1;
	vector <double> vectorOfSumElements(2 * polynomialDegree, 0.0);
	for (int k = 1; k <= 2 * polynomialDegree; k++)
		vectorOfSumElements[k - 1] = SumVectorElements(RaiseElementsOfVectorToPow(experimentalValuesOfX, k));
	vector <vector<double>> coefficientMatrix(polynomialDegree + 1, vector<double>(polynomialDegree + 1, 0.0));
	vector<double> vectorOfMultiplyParameters(polynomialDegree + 1, 0.0);
	for (int i = 0; i <= polynomialDegree; i++)
	{
		for (int j = 0; j <= polynomialDegree; j++)
		{
			if (i != 0 || j != 0)
			{
				int k = i + j - 1;
				coefficientMatrix[i][j] = vectorOfSumElements[k];
			}
			else
				coefficientMatrix[i][j] = numberOfExperimentalPoints;
		}
		vectorOfMultiplyParameters[i] = MultiplyParametersOfVectors(experimentalValuesOfY, RaiseElementsOfVectorToPow(experimentalValuesOfX, i));
	}
	vector<double> equationCoefficients = SolveGaussMethod(coefficientMatrix, vectorOfMultiplyParameters);
	return equationCoefficients;
}

double CalculateStandardDeviation(vector <double> experimentalValuesOfX, vector <double> experimentalValuesOfY, vector <double> equationCoefficients)
{
	double standardDeviation = 0;
	int numberOfExperimentalPoints = experimentalValuesOfX.size();
	int polynomialDegree = equationCoefficients.size() - 1;
	for (int i = 0; i < numberOfExperimentalPoints; i++)
	{
		double sumOfVectorElements = 0;
		for (int j = 0; j <= polynomialDegree; j++)
			sumOfVectorElements += equationCoefficients[j] * pow(experimentalValuesOfX[i], j);
		standardDeviation += pow(experimentalValuesOfY[i] - sumOfVectorElements, 2);
	}
	standardDeviation /= (numberOfExperimentalPoints - polynomialDegree - 1);
	return standardDeviation;
}

int main()
{
	vector <double> experimentalValuesOfX = { 2.40,3.50,5.00,6.89,10.00 }; // values of v in our equation: P = a+bv^2
	vector <double> experimentalValuesOfY = { 0.0141,0.0281,0.0562,0.1125,0.2250 }; // values of P in our equation
	vector <double> squaredValuesOfX = RaiseElementsOfVectorToPow(experimentalValuesOfX, 2); // we do equation P=a+bx, where x = v^2. reduction to a linear form of a function
	vector <double> equationCoefficients = SolveLSMApproximation(squaredValuesOfX, experimentalValuesOfY);
	cout << "equation: P=a+bv^2\nFound equation coefficients: \n";
	cout << "a = " << equationCoefficients[0] << ", b = " << equationCoefficients[1];
	double standardDeviation = CalculateStandardDeviation(experimentalValuesOfX, experimentalValuesOfY, equationCoefficients);
	cout << "\n\nStandard deviation of equations: \n" << standardDeviation;
	return 0;
}