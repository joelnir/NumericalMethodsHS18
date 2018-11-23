#include <cmath>
#include <iostream>
#include <functional>
#include <vector>

double GaussLegendre5(const std::function<double(double)> &f, double a, double b) {
	const std::vector<double> c = { -0.90617984593, -0.53846931010, 0.0, 0.53846931010, 0.90617984593 };
	const std::vector<double> w = { 0.23692688505, 0.47862867049, 0.56888888888, 0.47862867049, 0.23692688505 };

	// TODO: Task (a)
	// ...
	// ...
	// ...

	return .0; // dummy
}

double CompositeSimpson(const std::function<double(double)> &f, const std::vector<double> &x) {

	// TODO: Task (b)
	// ...
	// ...
	// ...

	return .0;
}

std::vector<double> LinSpace(int n, double a, double b) {
	std::vector<double> x(n);
	double d = (b - a) / (n - 1);
	for (int i = 0; i < n; ++i) {
		x[i] = i * d + a;
	}
	return x;
}

double f(double x) {
	return std::sqrt(x);
}

double F(double x) {
	double y = std::sqrt(x);
	return 2.0 / 3.0 * y * y * y;
}

int main() {
	int n = 5;
	double a = 0.0;
	double b = 1.0;

	std::cout << "Gauss-Legendre: " << GaussLegendre5(f, a, b) << std::endl;
	std::cout << "Simpson: " << CompositeSimpson(f, LinSpace(n, a, b)) << std::endl;
	std::cout << "Exact value: " << F(b) - F(a) << std::endl;

	return 0;
}

