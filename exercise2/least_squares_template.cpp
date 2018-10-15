#include <iostream>

#include <eigen3/Eigen/Dense>
#include <mgl2/mgl.h>

Eigen::MatrixXd Vandermonde(const Eigen::VectorXd &x, int n) {
    int m = x.size();
	Eigen::MatrixXd V(m, n);

    V.block(0, 0, m, 1) = Eigen::MatrixXd::Ones(m, 1);

    for(int i = 0; i < m; ++i){
        for(int j = 1; j < n; ++j){
            V(i, j) = V(i, j-1) * x(i);
        }
    }

	return V;
}

Eigen::VectorXd r(const Eigen::VectorXd &x) {
	return (1.0 / (1.0 + 25.0 * x.array() * x.array())).matrix();
}

int main() {
	int n = 11;				// Number of polynomial coefficients
	int m = 11;					// Number of samples
	Eigen::VectorXd x;		// Samples in [-1, 1]
    x.setLinSpaced(m, -1.0, 1.0);
	Eigen::MatrixXd V = Vandermonde(x, n);		// Vandermonde matrix
	Eigen::VectorXd y = r(x);		// r(x)
	Eigen::VectorXd a(n);	// Polynomial coefficients

    Eigen::VectorXd approx;

	Eigen::IOFormat PythonFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ";\n", "[", "]", "[", "]");

    std::cout << "Polynomial coefficients obtained by V" << std::endl;

	// Compute overfitted polynomial coefficients
	m = n;
	x.setLinSpaced(m, -1.0, 1.0);
    y = r(x);
    a = V.fullPivLu().solve(y);
	std::cout << "...overfitting:" << std::endl;
	std::cout << a.transpose().format(PythonFmt) << std::endl;

    approx = V*a;

    // Save data for plotting
    mglData x1_data, y1_data, x2_data, y2_data, r_data;
    double x1_arr[m], y1_arr[m], r_arr[m];
    for(int i = 0; i < m; ++i){
        x1_arr[i] = x(i);
        y1_arr[i] = approx(i);
        r_arr[i] = y(i);
    }
    x1_data.Set(x1_arr, m);
    y1_data.Set(y1_arr, m);
    r_data.Set(r_arr, m);


	// Compute least squares polynomial coefficients
	m = 3 * n;
	x.setLinSpaced(m, -1.0, 1.0);
    V = Vandermonde(x, n);
    y = r(x);
    a = V.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(y);

	std::cout << "...least squares:" << std::endl;
	std::cout << a.transpose().format(PythonFmt) << std::endl;

    approx = V*a;

    double x2_arr[m], y2_arr[m];
    for(int i = 0; i < m; ++i){
        x2_arr[i] = x(i);
        y2_arr[i] = approx(i);
    }
    x2_data.Set(x2_arr, m);
    y2_data.Set(y2_arr, m);

    mglGraph *gr = new mglGraph;
    gr->Title("Approximations");
  	gr->SetRanges(-1, 1, 0, 1);
  	gr->Axis();

  	gr->Plot(x1_data, y1_data,"r");
  	gr->AddLegend("Approximation 1 (inverse)","r");
  	gr->Plot(x2_data, y2_data,"g");
  	gr->AddLegend("Approximation 2 (LSQ)","g");
  	gr->Plot(x1_data, r_data,"b");
  	gr->AddLegend("Runge","b");

  	gr->Label('x',"X",0);
  	gr->Label('y', "Y",0);
    gr->Legend(3);
	gr->WriteFrame("matrix_comp.eps");
}

/*
 * d)
 *
 * In b) the solution is exactly the Runge function.
 * The least square solution from c) approximates the function,
 * showing an almost oscillating behaviour around the true value.
 *
 */
