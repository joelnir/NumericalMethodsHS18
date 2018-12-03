#include <cmath>
#include <iostream>
#include <functional>
#include <vector>

using namespace std;

double GaussLegendre5(const std::function<double(double)> &f, double a, double b) {
	const std::vector<double> c = { -0.90617984593, -0.53846931010, 0.0, 0.53846931010, 0.90617984593 };
	const std::vector<double> w = { 0.23692688505, 0.47862867049, 0.56888888888, 0.47862867049, 0.23692688505 };

    double sum = 0;

    function< double(double) > transform = [a, b](double x){
        return 0.5*(1-x)*a + 0.5*(x+1)*b;
    };

    for(int i = 0; i < 5; ++i){
        sum += w.at(i) * f(transform(c.at(i)));
    }

    return 0.5 * (b - a) * sum;
}

double CompositeSimpson(const std::function<double(double)> &f, const std::vector<double> &x) {

    int m = x.size() - 1;

    double sum = 0;

    sum += (1.0/6.0)*(x.at(1) - x.at(0))*f(x.at(0));

    for(int i = 1; i <= (m-1); ++i){
        sum += (1.0/6.0)*(x.at(i+1) - x.at(i-1))*f(x.at(i));
    }
    for(int i = 1; i <= m; ++i){
        sum += (2.0/3.0)*(x.at(i) - x.at(i-1))*f(0.5*(x.at(i) + x.at(i-1)));
    }

    sum += (1.0/6.0)*(x.at(m) - x.at(m-1))*f(x.at(m));

	return sum;
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

    cout << endl;

    cout << "5.2" << endl;

    auto exp_f = [](double x){return exp(x);};

	std::cout << "Gauss-Legendre: " << GaussLegendre5(exp_f, -3, 3) << std::endl;
	std::cout << "Exact value: " << (exp(3) - exp(-3)) << std::endl;


	return 0;
}

