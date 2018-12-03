#include <iostream>

#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

struct Newton {
	Newton(const Eigen::VectorXd &x) : _x(x), _a(x.size()) { }
	void Interpolate(const Eigen::VectorXd &y);
	double operator()(double x) const;

private:
	Eigen::VectorXd _x;	// nodes
	Eigen::VectorXd _a;	// coefficients
};

// Compute the coefficients in the Newton basis.
// O(n²), since matrix can be solved by abusing the fact that it is lower-triangular
void Newton::Interpolate(const Eigen::VectorXd &y) {
    int n = y.rows()-1;
    MatrixXd A = MatrixXd::Zero(n+1,n+1);

    // Set up system matrix A
    for(int row_i = 0; row_i <= n; ++row_i){
        for(int col_i = 0; col_i <= row_i; ++col_i){
            // 1 in first column
            if(col_i == 0){
                A(row_i,col_i) = 1;
            }
            else{
                A(row_i,col_i) = A(row_i,col_i-1)*(this->_x(row_i) - this->_x(col_i-1));
            }
        }
    }

    // Solve system

    this->_a = A.triangularView<Lower>().solve(y);
}

// Evaluate the interpolant at x.
double Newton::operator()(double x) const {
    double horner = 0;
    int n = this->_x.rows() - 1;

    horner = this->_a(n);
    for(int i = n-1; i >= 0; --i){
        horner = horner*(x - this->_x(i)) + this->_a(i);
    }

	return horner;
}

struct Lagrange {
	Lagrange(const Eigen::VectorXd &x);
	void Interpolate(const Eigen::VectorXd &y) { _y = y; }
	double operator()(double x) const;

private:
	Eigen::VectorXd _x;	// nodes
	Eigen::VectorXd _l;	// weights
	Eigen::VectorXd _y;	// coefficients
};

// Compute the weights l for given nodes x.
Lagrange::Lagrange(const Eigen::VectorXd &x) : _x(x), _l(x.size()), _y(x.size()) {
    int n = x.size() - 1;
    double prod;
    for(int i = 0; i <= n; ++i){
        prod = 1;
        for(int j = 0; j <= n; ++j){

            // For each except when i != j
            if(i != j){
                prod *= (x(i) - x(j));
            }

        }

        this->_l(i) = (1.0/prod);
    }
}

// Evaluate the interpolant at x.
// O(n), from calculation of w and final sum q
double Lagrange::operator()(double x) const {
    int n = this->_x.size() - 1;

    // Calculate w
    double w = 1;
    for(int i = 0; i <= n; ++i){
        w *= (x - this->_x(i));
    }

    // Calcualte Ls
    VectorXd L(n+1);
    for(int i = 0; i <= n; ++i){
        L(i) = w*(this->_l(i)/(x - this->_x(i)));
    }

    // Calculate final result
    double q = 0;
    for(int i = 0; i <= n; ++i){
        q += this->_y(i)*L(i);
    }

    return q;
}

// Runge function
Eigen::VectorXd r(const Eigen::VectorXd &x) {
	return (1.0 / (1.0 + 25.0 * x.array() * x.array())).matrix();
}

int main() {
	int n = 5;
	Eigen::VectorXd x;
	x.setLinSpaced(5, -1.0, 1.0);
	Eigen::VectorXd y = r(x);

	Newton p(x);
	p.Interpolate(y); // correct result: p._a = [0.0384615, 0.198939, 1.5252, -3.31565, 3.31565]

	Lagrange q(x);    // correct result: p._l = [0.666667, -2.66667, 4, -2.66667, 0.666667]
	q.Interpolate(y);

	// Compute difference of p and q.
	int m = 22;
	double offset = 0.08333333333;
	x.setLinSpaced(m, -1.0 + offset, 1.0 - offset);
	double norm2 = .0;
	for (int i = 0; i < m; ++i) {
		double d = p(x(i)) - q(x(i));
		norm2 += d * d;
	}

	// By uniquenss of the interpolation polynomial, we expect p = q.
	std::cout << "This number should be close to zero: " << norm2 << std::endl;

	return 0;
}

