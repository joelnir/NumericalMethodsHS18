#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

using namespace std;

std::vector<double> Heun(const std::function<double(double, double)> &f, const std::vector<double> &t, double y0) {
	int n = t.size();
	std::vector<double> y(n);

    // Initial value
    y[0] = y0;

    double h = t[1] - t[0];

    double k1, k2;
    for(int k = 1; k < n; ++k){
        k1 = f(t[k-1], y[k-1]);
        k2 = f(t[k], y[k-1] + h*k1);
        y[k] = y[k-1] + (h/2)*(k1 + k2);
    }

    return y;
}

// right-hand side
double f(double t, double y) {
	return std::exp(-2. * t) - 2. * y;
}

// Exact solution for Y(0) = 0.
double Y(double t) {
	return t * std::exp(-2. * t);
}

// Create an evenly spaced grid of n points on the interval [a,b].
std::vector<double> LinSpace(int n, double a, double b) {
	std::vector<double> x(n);
	double d = (b - a) / (n - 1);
	for (int i = 0; i < n; ++i) {
		x[i] = i * d + a;
	}
	return x;
}

template <class T>
std::ostream & operator<<(std::ostream &os, const std::vector<T> &v) {
	os << "[";
	if (!v.empty()) {
		os << v[0];
		for (int i = 1; i < v.size(); ++i) os << ", " << v[i];
	}
    os << "]";
    return os;
}

int main() {
    double T = 1.;

    // A simple test:
	std::cout << "Heun:  " << Heun(f, LinSpace(16, .0, T), .0).back() << std::endl;
	std::cout << "exact: " << Y(T) << std::endl;

	// Convergence analysis:
    std::vector<int> n = { 4, 8, 16, 32 };
    int J = n.size();
    std::vector<double> h(J);
    std::vector<double> error(J);
    for (int j = 0; j < J; ++j) {
    	h[j] = 1. / n[j];
    	std::vector<double> t = LinSpace(n[j] + 1, .0, T);
    	std::vector<double> y = Heun(f, t, .0);

        error[j] = (abs(Y(1.0) - y[n[j]]));
    }

    std::cout << "h = " << h << std::endl;
    std::cout << "error = " << error << std::endl;

    return 0;
}
