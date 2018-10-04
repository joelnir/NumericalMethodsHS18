// compile with: g++ template_structured_matrix_vector.cpp -I/usr/include/eigen3 -lmgl

//#include <chrono>
#include <iostream>
#include <iomanip>
//#include <limits>
//#include <ratio>
#include <vector>
#include <cmath>
#include <time.h>

#include <Eigen/Dense>
#include <mgl2/mgl.h>

using namespace Eigen;

/* \brief compute $\mathbf{A}\mathbf{x}$
 * \mathbf{A} is defined by $(\mathbf{A})_{i,j} := \min {i,j}$
 * \param[in] x vector x for computation of A*x = y
 * \param[out] y = A*x
 */
void multAminSlow(const VectorXd & x, VectorXd & y) {
    unsigned int n = x.size();

    VectorXd one = VectorXd::Ones(n);
    VectorXd linsp = VectorXd::LinSpaced(n,1,n);
    y = ( ( one * linsp.transpose() )
          .cwiseMin( linsp * one.transpose()) ) * x;
}
/*
 * a)
 *
 * The time complexity of the code is in O(n^2).
 */

/* \brief compute $\mathbf{A}\mathbf{x}$
 * \mathbf{A} is defined by $(\mathbf{A})_{i,j} := \min {i,j}$
 * Instead of a "Matlab style" construcion of the product,
 * we use simple loops.
 * \param[in] x vector x for computation of A*x = y
 * \param[out] y = A*x
 */
void multAminLoops(const VectorXd & x, VectorXd & y) {
    unsigned int n = x.size();

    MatrixXd A(n,n);

    for(unsigned int i = 0; i < n; ++i) {
        for(unsigned int j = 0; j < n; ++j) {
            A(i,j) = std::min(i+1,j+1);
        }
    }
    y = A * x;
}

/* \brief compute $\mathbf{A}\mathbf{x}$
 * This function has optimal complexity.
 * \mathbf{A} is defined by $(\mathbf{A})_{i,j} := \min {i,j}$
 * \param[in] x vector x for computation of A*x = y
 * \param[out] y = A*x
 */
void multAmin(const VectorXd & x, VectorXd & y) {
    unsigned int n = x.size();
    VectorXd part_sums = VectorXd(n);

    double sum = 0;
    //Calculate partial sums
    for(int i = (n-1); i >= 0; --i){
        sum += x(i);
        part_sums(i) = sum;
    }

    y.resize(n);

    if(n > 0){
        y(0) = part_sums(0);
    }

    for(int i = 1; i < n; ++i){
        y(i) = y(i - 1) + part_sums(i);
    }
}

int main(void) {
    // Testing correctness of the code
    unsigned int M = 10;
    VectorXd xa = VectorXd::Random(M);
    VectorXd ys, yf;

    multAmin(xa, yf);
    multAminSlow(xa, ys);
    // Error should be small
    std::cout << "||ys-yf|| = " << (ys - yf).norm() << std::endl;


    unsigned int nLevels = 9;
	unsigned int *n = new unsigned int[nLevels];
	double *minTime = new double[nLevels];
	double *minTimeLoops = new double[nLevels];
	double *minTimeEff = new double[nLevels];

	n[0] = 4;
	for (unsigned int i=1; i<nLevels; i++){
		n[i] = 2*n[i-1];
    }

    /*
     * c)
     *
     * Runtimes (s):
     *
     *   n  | slow    | fast
     * -----|---------|---------
     * 2^4  |0.580e-4 |0.600e-5
     * 2^5  |0.154e-3 |0.130e-4
     * 2^6  |0.529e-3 |0.240e-4
     * 2^7  |0.200e-2 |0.480e-4
     * 2^8  |0.686e-2 |0.910e-4
     * 2^9  |0.233e-1 |0.144e-3
     * 2^10 |0.726e-1 |0.252e-3
     *
     */
    clock_t start_time, slow_time, fast_time;
    VectorXd y_time;

    std::cout.precision(10);
    std::cout << std::endl;

    bool run_timings = false;
    if(run_timings){
        for(int i = 4; i <= 10; ++i){
            int n = std::pow(2, i);

            std::cout << "Runtime for n = 2 ^ " << i << std::endl;
            std::cout << "---------------------" << std::endl;

            VectorXd test_vec = VectorXd::Random(n);

            // Time slow function
            start_time = clock();
            multAminSlow(test_vec, y_time);
            slow_time = clock() - start_time;

            // Time fast function
            start_time = clock();
            multAmin(test_vec, y_time);
            fast_time = clock() - start_time;

            std::cout << "Slow: " << ((float)slow_time)/CLOCKS_PER_SEC << " s" << std::endl;
            std::cout << "Fast: " << ((float)fast_time)/CLOCKS_PER_SEC << " s" << std::endl;
            std::cout << std::endl;
        }
    }


    // Plotting with MathGL
    double nMgl[nLevels];
    double ref1[nLevels], ref2[nLevels];
    for (int i=0; i<nLevels; i++) {
    	nMgl[i] = n[i];
    	ref1[i] = 1e-8*pow(n[i],2);
    	ref2[i] = 1e-7*n[i];
    }

    mglData matSize;
    matSize.Link(nMgl, nLevels);

    mglData data1, data2;
    mglData dataRef1, dataRef2;
  	data1.Link(minTime, nLevels);
  	data2.Link(minTimeEff, nLevels);
  	dataRef1.Link(ref1,nLevels);
  	dataRef2.Link(ref2,nLevels);

  	mglGraph *gr = new mglGraph;
    gr->Title("Runtime of multAmin");
  	gr->SetRanges(n[0],n[0]*pow(2,nLevels-1),1e-6,1e+1);  gr->SetFunc("lg(x)","lg(y)");
  	gr->Axis();
  	gr->Plot(matSize,data1,"k +"); gr->AddLegend("slow","k +");
  	gr->Plot(matSize,data2,"r +"); gr->AddLegend("efficient","r +");
  	gr->Plot(matSize,dataRef1,"k"); gr->AddLegend("O(n^2)","k");
  	gr->Plot(matSize,dataRef2,"r"); gr->AddLegend("O(n)","r");
  	gr->Label('x',"Matrix size [n]",0);
  	gr->Label('y', "Runtime [s]",0);
    gr->Legend(2);
	gr->WriteFrame("multAmin_comparison.eps");


    // The following code is just for demonstration purposes.
    // Build Matrix B with dimension 10x10
    unsigned int nn = 10;
    MatrixXd B = MatrixXd::Zero(nn,nn);
    for(unsigned int i = 0; i < nn; ++i) {
        B(i,i) = 2;
        if(i < nn-1) B(i+1,i) = -1;
        if(i > 0) B(i-1,i) = -1;
    }
    B(nn-1,nn-1) = 1;

    /*
     * d)
     *
     * (2  -1 0  0  -  0)
     * (-1 2  -1 0  -  0)
     * (0  -1 2  -1 -  0)
     * (0  0  -1 2  -  0)
     * (|  |  |  |  |  |)
     * (0  0  0  0  -  1)
     */

    /*
     * e)
     *
     * The computation gives e_j as a result, showing that B = A^-1 (inverse).
     *
     */
    int j = 4;
    VectorXd e_j = VectorXd::Zero(nn);
    e_j(j) = 1;

    VectorXd res = VectorXd();
    multAmin((B*e_j), res);

    std::cout << "Result of ABe_j: " << std::endl << res << std::endl;
}

