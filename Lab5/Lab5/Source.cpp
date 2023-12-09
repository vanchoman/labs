#include <iostream>
#include <cmath>
#include <vector>

double calculateFunction2(const double& x, const double& y) {
    return (std::pow(x, 2) + 2 * y);
}

double calculateIntegralForCube(const std::vector<double>& x, const std::vector<double>& y, int n, int m, double step1, double step2) {
    double integral = 0.0;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            integral += calculateFunction2(x[2 * i], y[2 * j]) +
                calculateFunction2(x[2 * i + 2], y[2 * j]) +
                calculateFunction2(x[2 * i + 2], y[2 * j + 2]) +
                calculateFunction2(x[2 * i], y[2 * j + 2]) +
                4 * (calculateFunction2(x[2 * i + 1], y[2 * j]) +
                    calculateFunction2(x[2 * i + 2], y[2 * j + 1]) +
                    calculateFunction2(x[2 * i + 1], y[2 * j + 2]) +
                    calculateFunction2(x[2 * i], y[2 * j + 1])) +
                16 * calculateFunction2(x[2 * i + 1], y[2 * j + 1]);
        }
    }
    return integral * step1 * step2 / 9;
}

double calculateCubeSimpsonIntegral(const std::vector<double>& firstSpan, const std::vector<double>& secondSpan, double eps) {
    int n = 4, m = n;
    double step1 = (firstSpan[1] - firstSpan[0]) / (2 * n);
    double step2 = (secondSpan[1] - secondSpan[0]) / (2 * m);

    std::vector<double> x(2 * (n + 2)), y(2 * (m + 2));
    for (int i = 0; i <= 2 * (n + 1); ++i) {
        x[i] = firstSpan[0] + i * step1;
    }
    for (int i = 0; i <= 2 * (m + 1); ++i) {
        y[i] = secondSpan[0] + i * step2;
    }

    return calculateIntegralForCube(x, y, n, m, step1, step2);
}

int main() {
    double eps = 1e-5;

    std::vector<double> firstSpan = { 0, 2.0 };
    std::vector<double> secondSpan = { 0, 1.0 };

    double cubeIntegral = calculateCubeSimpsonIntegral(firstSpan, secondSpan, eps);
    std::cout << "\nCube integral by Simpson: " << cubeIntegral << std::endl;

    return 0;
}