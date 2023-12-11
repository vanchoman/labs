#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

void OutMatrixInConsole(const vector<vector<double>>& A, const vector<double>& b);

enum UserChoiceAboutSolMethod
{
	GaussMethod = 1,
	LDLFactorization = 2,
};

enum UserChoiceAboutInputMethod
{
	manualEntryMatrix = 1,
	defaultMatrix = 2,
};

bool CheckRowProportionality(const vector<double>& firstRows, const vector<double>& secondRows) // checking rows for proportionality
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

vector<double> MultiplyMatrixVector(vector<vector<double>> matrixOfElements, vector<double> vectorOfElements)
{
	vector<double> resultVector(matrixOfElements.size());
	int vectorSize = resultVector.size();
	for (int i = 0; i < vectorSize; i++)
	{
		resultVector[i] = 0;
		for (int j = 0; j < vectorSize; j++)
			resultVector[i] += matrixOfElements[i][j] * vectorOfElements[j];
	}
	return resultVector;
}

bool CheckProportionality(const vector<vector<double>>& A, const vector<double>& b) // checking a matrix for row proportionality
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

double FindMaxInRV(const vector<vector<double>>& A, const vector<double>& b, const vector<double>& firstRoots) // finding the residual vector and its norm
{
	vector<double> vectorCalculationDiff(firstRoots.size());
	double maxInVCD = 0;
	int vectorSize = firstRoots.size();
	for (int i = 0; i < vectorSize; i++)
	{
		for (int j = 0; j < vectorSize; j++)
			vectorCalculationDiff[i] += A[i][j] * firstRoots[j];
		vectorCalculationDiff[i] -= b[i];
	}
	maxInVCD = abs(vectorCalculationDiff[0]);
	for (int i = 1; i < vectorSize; i++)
	{
		if (maxInVCD < abs(vectorCalculationDiff[i]))
			maxInVCD = abs(vectorCalculationDiff[i]);
	}
	return  maxInVCD;
}

double CalculateError(const vector<double>& firstRoots, const vector<double>& secondRoots) // calculation of relative error
{
	int vectorSize = firstRoots.size();
	double calcError = 0, maxDiff = 0, maxRoots = 0;
	for (int i = 0; i < vectorSize; i++)
	{
		if (abs(secondRoots[i]) - abs(firstRoots[i]) > maxDiff)
			maxDiff = abs(secondRoots[i]) - abs(firstRoots[i]);
		if (firstRoots[i] > maxRoots)
			maxRoots = firstRoots[i];
	}
	calcError = abs(maxDiff) / abs(maxRoots);
	return calcError;
}

void OutMatrixInConsole(const vector<vector<double>>& A, const vector<double>& b) // output matrix to console
{
	int vectorSize = A.size();
	for (int i = 0; i < vectorSize; i++)
	{
		for (int j = 0; j < vectorSize; j++)
			cout << A[i][j] << setw(11);
		cout << setw(3) << "| " << b[i] << endl;
	}
}

vector<double> SolveSystemWithLDLFactorization(vector<vector<double>>& A, vector<double>& b)
{
	int sizeOfMatrix = A.size();
	vector<double> diagonalMatrix(sizeOfMatrix), vectorOfRoots(sizeOfMatrix), y(sizeOfMatrix), z(sizeOfMatrix);
	vector<vector<double>> lowerTriangularMatrix(sizeOfMatrix, vector<double>(sizeOfMatrix, 0));
	vector<vector<double>> upperTriangularMatrix(lowerTriangularMatrix);
	for (int i = 0; i < sizeOfMatrix; i++) // make a unit diagonal in matrices
		lowerTriangularMatrix[i][i] = upperTriangularMatrix[i][i] = 1;
	for (int j = 0; j < sizeOfMatrix; j++)
	{
		double sumofMultElements = 0;
		for (int i = 0; i < j; i++)
			sumofMultElements += diagonalMatrix[i] * pow(lowerTriangularMatrix[j][i], 2);
		diagonalMatrix[j] = A[j][j] - sumofMultElements; // filling out our diagonal matrix
		for (int i = j + 1; i < sizeOfMatrix; i++)
			lowerTriangularMatrix[i][j] = (A[i][j] - sumofMultElements) / diagonalMatrix[j]; // fill the lower diagonal matrix
	}
	for (int i = 0; i < sizeOfMatrix - 1; i++)
	{
		for (int j = i + 1; j < sizeOfMatrix; j++)
			upperTriangularMatrix[i][j] = lowerTriangularMatrix[j][i];
	}
	for (int i = 0; i < sizeOfMatrix; i++)
	{
		double sumofMultElements = 0;
		for (int j = 0; j < i; j++)
			sumofMultElements += lowerTriangularMatrix[i][j] * y[j];
		y[i] = b[i] - sumofMultElements;
		z[i] = y[i] / diagonalMatrix[i];
	}
	for (int i = sizeOfMatrix - 1; i >= 0; i--)
	{
		double sumofMultElements = 0;
		for (int j = sizeOfMatrix - 1; j > i; j--)
			sumofMultElements += upperTriangularMatrix[i][j] * vectorOfRoots[j];
		vectorOfRoots[i] = z[i] - sumofMultElements;
	}
	return vectorOfRoots;
}

int main()
{
	vector<vector<double>> A;
	vector<double> b;
	int userChoice;
	cout << "what do you want to check?\n" <<
		(int)UserChoiceAboutSolMethod::GaussMethod << " - Gaussian method\n" <<
		(int)UserChoiceAboutSolMethod::LDLFactorization << " - LDL factorization\n";
	cin >> userChoice;
	switch ((UserChoiceAboutSolMethod)userChoice)
	{
	case UserChoiceAboutSolMethod::GaussMethod:
	{
		int userChoice;
		cout << "\nhow do you want to initialize the matrix?\n"
			<< (int)UserChoiceAboutInputMethod::manualEntryMatrix << " - enter the matrix manually\n"
			<< (int)UserChoiceAboutInputMethod::defaultMatrix << " - use an already defined matrix\n"
			<< "any other value terminates the program!\n\nyour answer: ";
		cin >> userChoice;
		switch ((UserChoiceAboutInputMethod)userChoice)
		{
		case UserChoiceAboutInputMethod::manualEntryMatrix:
		{
			cout << "enter matrix order: ";
			int matrixSize;
			cin >> matrixSize;
			A.assign(matrixSize, vector<double>(matrixSize));
			b.resize(matrixSize);
			for (int i = 0; i < matrixSize; i++)
			{
				for (int j = 0; j < matrixSize; j++)
				{
					cout << "enter element with index A[" << i + 1 << "][" << j + 1 << "]:\n";
					cin >> A[i][j];
				}
			}
			for (int i = 0; i < matrixSize; i++)
			{
				cout << "enter element b[" << i + 1 << "]:\n";
				cin >> b[i];
			}
			break;
		}

		case UserChoiceAboutInputMethod::defaultMatrix:
		{
			A = { {0.14,0.24,-0.84},
				{1.07,-0.83,0.56},
				{0.64,0.43,-0.38} };
			b = { 1.11,0.48,-0.83 };
			break;
		}

		default:
			cout << "the program has ended\n";
			return 0;
		}

		cout << "\nentered matrix: \n";
		OutMatrixInConsole(A, b);
		if (CheckProportionality(A, b))
		{
			cout << "\nrows are proportional, the roots of the system cannot be calculated\n";
			return 0;
		}
		cout << "\nanswer: \n";
		vector<double> firstRoots = SolveGaussMethod(A, b);
		for (int i = 0; i < firstRoots.size(); i++)
			cout << " x" << i + 1 << " = " << firstRoots[i];
		double maxNormaOfRV = FindMaxInRV(A, b, firstRoots);
		cout << "\n\nnorma of residual vector: \n" << maxNormaOfRV << endl;
		vector<double> secondRoots = SolveGaussMethod(A, MultiplyMatrixVector(A, firstRoots));
		cout << "\nsecond solution: \n";
		for (int i = 0; i < secondRoots.size(); i++)
			cout << " x" << i + 1 << " = " << secondRoots[i];
		cout << "\n\ncalculation error: \n";
		double calculationError = CalculateError(firstRoots, secondRoots);
		cout << calculationError << endl;
		break;
	}

	case UserChoiceAboutSolMethod::LDLFactorization:
	{
		int userChoice1;
		vector<double> vectorOfEigenvalues(3);
		cout << "you need to enter three eigenvalues or or you can use default eigenvalues\n"
			<< "1 - use default\n2 - enter manually\nyour answer: ";
		cin >> userChoice1;
		if (userChoice1 == 1)
		{
			vectorOfEigenvalues[0] = 1;
			vectorOfEigenvalues[1] = pow(10, 3);
			vectorOfEigenvalues[2] = pow(10, 6);
		}
		else if (userChoice1 == 2)
		{
			for (int i = 0; i < vectorOfEigenvalues.size(); i++)
			{
				cout << i + 1 << " lambda: ";
				cin >> vectorOfEigenvalues[i];
			}
		}
		else
		{
			cout << "Enter correct choice!\n";
			main();
		}
		vector<vector<double>> A(3, vector<double>(3));
		A[0][0] = 2 * vectorOfEigenvalues[0] + 4 * vectorOfEigenvalues[1];
		A[1][1] = 2 * vectorOfEigenvalues[0] + vectorOfEigenvalues[1] + 3 * vectorOfEigenvalues[2];
		A[2][2] = 2 * vectorOfEigenvalues[0] + vectorOfEigenvalues[1] + 3 * vectorOfEigenvalues[2];
		A[0][1] = A[1][0] = 2 * (vectorOfEigenvalues[0] - vectorOfEigenvalues[1]);
		A[2][0] = A[0][2] = 2 * (vectorOfEigenvalues[0] - vectorOfEigenvalues[1]);
		A[1][2] = A[2][1] = 2 * vectorOfEigenvalues[0] + vectorOfEigenvalues[1] - 3 * vectorOfEigenvalues[2];
		vector<double> b = { -4 * vectorOfEigenvalues[0] - 2 * vectorOfEigenvalues[1], -4 * vectorOfEigenvalues[0] + vectorOfEigenvalues[1] + 9 * vectorOfEigenvalues[2], -4 * vectorOfEigenvalues[0] + vectorOfEigenvalues[1] - 9 * vectorOfEigenvalues[2] };
		OutMatrixInConsole(A, b);
		vector<double> vectorOfRoots = SolveSystemWithLDLFactorization(A, b);
		cout << "roots found using LDL factorization: " << endl;
		int vectorSize = vectorOfRoots.size();
		for (int i = 0; i < vectorSize; i++)
			cout << " x" << i + 1 << " = " << vectorOfRoots[i];
	}
	}
	return 0;
}
