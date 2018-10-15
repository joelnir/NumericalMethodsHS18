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
	int m;					// Number of samples
	Eigen::VectorXd x;		// Samples in [-1, 1]
	Eigen::MatrixXd V;		// Vandermonde matrix
	Eigen::VectorXd y;		// r(x)
	Eigen::VectorXd a1(n);	// Polynomial coefficients
	Eigen::VectorXd a2(n);	// Polynomial coefficients

	Eigen::IOFormat PythonFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ";\n", "[", "]", "[", "]");

    std::cout << "Polynomial coefficients obtained by V" << std::endl;

	// Compute overfitted polynomial coefficients
	m = n;
	x.setLinSpaced(m, -1.0, 1.0);
    y = r(x);
	V = Vandermonde(x, n);
    a1 = V.fullPivLu().solve(y);
	std::cout << "...overfitting:" << std::endl;
	std::cout << a1.transpose().format(PythonFmt) << std::endl;

	// Compute least squares polynomial coefficients
	m = 3 * n;
	x.setLinSpaced(m, -1.0, 1.0);
    V = Vandermonde(x, n);
    y = r(x);
    a2 = V.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(y);

	std::cout << "...least squares:" << std::endl;
	std::cout << a2.transpose().format(PythonFmt) << std::endl;

    // Save data for plotting
    int p = 30;
    Eigen::VectorXd x_p(p);
    x_p.setLinSpaced(p, -1.0, 1.0);

    Eigen::VectorXd y1 = Vandermonde(x_p, n)*a1;
    Eigen::VectorXd y2 = Vandermonde(x_p, n)*a2;
    Eigen::VectorXd yr = r(x_p);

    mglData x_data, y1_data, y2_data, r_data;
    double x_arr[p], y1_arr[p], r_arr[p], y2_arr[p];
    for(int i = 0; i < p; ++i){
        x_arr[i] = x_p(i);
        y1_arr[i] = y1(i);
        y2_arr[i] = y2(i);
        r_arr[i] = yr(i);
    }
    x_data.Set(x_arr, m);
    y1_data.Set(y1_arr, m);
    r_data.Set(r_arr, m);
    y2_data.Set(y2_arr, m);

    mglGraph *gr = new mglGraph;
    gr->Title("Approximations");
  	gr->SetRanges(-1.0, 1.0, 0.0, 1.0);
  	gr->Axis();

  	gr->Plot(x_data, y1_data,"r");
  	gr->AddLegend("Approximation 1 (inverse)","r");
  	gr->Plot(x_data, y2_data,"g");
  	gr->AddLegend("Approximation 2 (LSQ)","g");
  	gr->Plot(x_data, r_data,"b");
  	gr->AddLegend("Runge","b");

  	gr->Label('x',"X",0);
  	gr->Label('y', "Y",0);
    gr->Legend(3);
	gr->WriteFrame("matrix_comp.eps");
}

/*
 * d)
 *
 * In b) the solution has overfit to the data points from the runge function.
 * The least square solution from c) approximates the function,
 * showing an almost oscillating behaviour around the true value.
 *
 */
