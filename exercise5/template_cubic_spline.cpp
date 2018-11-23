// to compile run: g++ -std=gnu++11 cubic_spline.cpp -lmgl
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <mgl2/mgl.h>

using namespace Eigen;

MatrixXd cubicSpline(const VectorXd &T, const VectorXd &Y) {
	// returns the matrix representing the spline interpolating the data
	// with abscissae T and ordinatae Y. Each column represents the coefficients
	// of the cubic polynomial on a subinterval.
	// Assumes T is sorted, has no repeated elements and T.size() == Y.size().

	int n = T.size() - 1; // T and Y have length n+1

	// TODO: build the spline matrix with polynomials' coefficients
	MatrixXd spline(4, n);

	return spline;
}

VectorXd evalCubicSpline(const MatrixXd &S, const VectorXd &T, const VectorXd &evalT) {
	// Returns the values of the spline S calculated in the points X.
	// Assumes T is sorted, with no repetetions.

	int n = evalT.size();
	VectorXd out(n);

	// TODO: fill out

	return out;
}

int main() {
	// tests
	VectorXd T(9);
	VectorXd Y(9);
	T << 0, 0.4802, 0.7634, 1, 1.232, 1.407, 1.585, 1.879, 2;
	Y << 0., 0.338, 0.7456, 0, -1.234, 0 , 1.62, -2.123, 0;

	int len = 1 << 9;
	VectorXd evalT = VectorXd::LinSpaced(len, T(0), T(T.size()-1));

	VectorXd evalSpline = evalCubicSpline(cubicSpline(T, Y), T, evalT);

 	mglData datx, daty;
	datx.Link(evalT.data(), len);
	daty.Link(evalSpline.data(), len);
	mglGraph gr;
	gr.SetRanges(0, 2, -3, 3);
	gr.Plot(datx, daty, "0");
	gr.WriteFrame("spline.eps");
}
	
