#include<bits/stdc++.h>
#include<eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

MatrixXf mult(const MatrixXf &A, const MatrixXf &B) {
	int N = A.rows();
	MatrixXf C = MatrixXf::Zero(N, N);

	//TODO: Point (a)

	return C;
}

MatrixXf mult_rec(const MatrixXf &A, const MatrixXf &B) {
	int N = A.rows();
	MatrixXf C(N, N);

	//TODO: Point (c)

	return C;
}
	

MatrixXf strassen(const MatrixXf &A, const MatrixXf &B) {
	int N = A.rows();
	MatrixXf C(N, N);

	//TODO: Point (e)

	return C;
}

int main() {
	srand(time(0));
	cout << setprecision(6) << setfill(' ');

	for (int i = 1; i < 9; i++) {
		int N = 1 << i;		
		cout << "Matrix size = " << N << endl;
		MatrixXd AA = MatrixXd::Random(N, N);
		MatrixXd BB = MatrixXd::Random(N, N);
		MatrixXd ans = AA*BB;
		MatrixXf A = AA.cast<float>();
		MatrixXf B = BB.cast<float>();

		auto start = std::chrono::steady_clock::now();
		MatrixXf W = mult(A, B);
		auto finish = std::chrono::steady_clock::now();
		cout << setw(24) << " " <<  setw(15) 
			<< "Time (s)" << setw(20) << "Error (l2-norm)"  << endl;
		cout << setw(24) << "Naive iterative "<< setw(15) 
			<< std::chrono::duration_cast<std::chrono::duration<double> > 
			(finish - start).count() << setw(20) << (W.cast<double>() - ans).norm() << endl;
	
		start = std::chrono::steady_clock::now();
		MatrixXf X = mult_rec(A, B);
		finish = std::chrono::steady_clock::now();
		cout << setw(24) << "Naive recursive "<< setw(15) 
			<< std::chrono::duration_cast<std::chrono::duration<double> > 
			(finish - start).count() << setw(20) << (X.cast<double>() - ans).norm() << endl;
		start = std::chrono::steady_clock::now();
		MatrixXf Y = strassen(A, B);
		finish = std::chrono::steady_clock::now();
		cout << setw(24) << "Strassen recursive " << setw(15) 
			<< std::chrono::duration_cast<std::chrono::duration<double> >
			(finish - start).count() << setw(20) << (Y.cast<double>() - ans).norm() << endl;
		start = std::chrono::steady_clock::now();	
		MatrixXf Z = A*B;
		finish = std::chrono::steady_clock::now();
		cout << setw(24) << "Eigen built-in "<< setw(15) 
			<< std::chrono::duration_cast<std::chrono::duration<double> >
			(finish - start).count() << setw(20) << (Z.cast<double>() - ans).norm() << "\n\n\n";
	}
}

