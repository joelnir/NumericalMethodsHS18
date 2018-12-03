// to compile run: g++ -std=gnu++11 cubic_spline.cpp -lmgl
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <mgl2/mgl.h>

using namespace Eigen;
using namespace std;

MatrixXd cubicSpline(const VectorXd &T, const VectorXd &Y) {
	// returns the matrix representing the spline interpolating the data
	// with abscissae T and ordinatae Y. Each column represents the coefficients
	// of the cubic polynomial on a subinterval.
	// Assumes T is sorted, has no repeated elements and T.size() == Y.size().

	int n = T.size() - 1; // T and Y have length n+1
    VectorXd h = T.bottomRows(n) - T.topRows(n);

    // Construct r vector
    VectorXd r(n-1);
    for(int i = 1; i <= n-1; ++i){
        r(i-1) = ((Y(i+1) - Y(i))/h(i)) -((Y(i) - Y(i-1))/h(i-1));
    }

    MatrixXd system(n-1, n-1);

    // Set up system matrix
    system(0,0) = (h(0) + h(1))/3;
    system(0,2) = h(1)/6;

    for(int row = 1; row <= n-3; ++row){
        system(row, row-1) = h(row)/6;
        system(row, row) = (h(row) + h(row+1))/3;
        system(row, row+1) = h(row+1)/6;
    }

    system(n-2,n-3) = h(n-2)/6;
    system(n-2,n-2) = (h(n-2) + h(n-1))/3;

    // Solve system to get sigmas
    VectorXd sigma = VectorXd::Zero(n+1);

    sigma.block(1,0,n-1,1) = system.fullPivLu().solve(r);

    // Fill spline array with coefficients
	MatrixXd spline(4, n);
    for(int i = 0; i < n; ++i){
        spline(0, i) = Y(i); //a
        spline(1, i) = ((Y(i+1) - Y(i))/h(i)) - ((h(i)*(2*sigma(i) + sigma(i+1) ))/6);//b
        spline(2, i) = sigma(i)/2; //c
        spline(3, i) = (sigma(i+1) - sigma(i))/(6*h(i)); //d
    }

	return spline;
}

VectorXd evalCubicSpline(const MatrixXd &S, const VectorXd &T, const VectorXd &evalT) {
	// Returns the values of the spline S calculated in the points T.
	// Assumes T is sorted, with no repetetions.

	int n = evalT.size();
	VectorXd out(n);

    double evaluation;
    VectorXd v(4);
    VectorXd coeff(4);
    for(int eval_i = 0; eval_i < n; ++eval_i){

        // Find interval to evaluate in
        int interval_i = 0;
        while(evalT(eval_i) > T(interval_i + 1)){
            interval_i++;
        }

        // Create vector of 1,t,t²,t³
        v(0) = 1;
        for(int i = 1; i <= 3; ++i){
            v(i) = (evalT(eval_i) - T(interval_i))*v(i-1);
        }

        // Evaluate using coefficients from spline
        coeff = S.block(0,interval_i,4,1);
        evaluation = coeff.dot(v);
        out(eval_i) = evaluation;
    }

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

    MatrixXd splineC = cubicSpline(T, Y);
	VectorXd evalSpline = evalCubicSpline(splineC, T, evalT);

 	mglData datx, daty;
	datx.Link(evalT.data(), len);
	daty.Link(evalSpline.data(), len);
	mglGraph gr;
	gr.SetRanges(0, 2, -3, 3);
	gr.Plot(datx, daty, "0");
	gr.WriteFrame("spline.eps");
}

