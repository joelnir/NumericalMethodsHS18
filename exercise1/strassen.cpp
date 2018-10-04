#include<bits/stdc++.h>
#include<eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

/*
 * a)
 */
MatrixXf mult(const MatrixXf &A, const MatrixXf &B) {
	int N = A.rows();
	MatrixXf C = MatrixXf::Zero(N, N);

    float new_c = 0;

    for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
            new_c = 0;

            for(int k = 0; k < N; ++k){
                new_c += A(i,k)*B(k, j);
            }

            C(i, j) = new_c;
        }
    }

	return C;
}

/*
 * b) O(N^3)
 */

/*
 * c)
 */
MatrixXf mult_rec(const MatrixXf &A, const MatrixXf &B) {
	int N = A.rows();
	MatrixXf C(N, N);

    if(N == 1){
        // Base case A, B = 1x1 matrix
        C(0,0) = A(0,0)*B(0,0);
    }
    else{
        int block_size = N/2;
        for(int i = 0; i < 2; ++i){
            for(int j = 0; j < 2; ++j){
                C.block(i*block_size, j*block_size, block_size, block_size) =
                    mult_rec(
                            A.block(i*block_size, 0, block_size, block_size),
                            B.block(0, j*block_size, block_size, block_size)
                            ) +
                    mult_rec(
                            A.block(i*block_size, block_size, block_size, block_size),
                            B.block(block_size, j*block_size, block_size, block_size)
                            );
            }
        }
    }

	return C;
}

/*
 * d)
 *
 * O(N^3)
 *
 */

/*
 * e)
 */
MatrixXf strassen(const MatrixXf &A, const MatrixXf &B) {
	int N = A.rows();
	MatrixXf C(N, N);

    if(N == 1){
        // Base case A, B = 1x1 matrix
        C(0,0) = A(0,0)*B(0,0);
    }
    else{
        int block_size = N/2;

        MatrixXf A_11 = A.block(0, 0, block_size, block_size);
        MatrixXf A_12 = A.block(0, block_size, block_size, block_size);
        MatrixXf A_21 = A.block(block_size, 0, block_size, block_size);
        MatrixXf A_22 = A.block(block_size, block_size, block_size, block_size);

        MatrixXf B_11 = B.block(0, 0, block_size, block_size);
        MatrixXf B_12 = B.block(0, block_size, block_size, block_size);
        MatrixXf B_21 = B.block(block_size, 0, block_size, block_size);
        MatrixXf B_22 = B.block(block_size, block_size, block_size, block_size);

        MatrixXf M_1 = strassen((A_11+A_22), (B_11+B_22));
        MatrixXf M_2 = strassen((A_21+A_22), B_11);
        MatrixXf M_3 = strassen(A_11, (B_12-B_22));
        MatrixXf M_4 = strassen(A_22, (B_21-B_11));
        MatrixXf M_5 = strassen((A_11+A_12), (B_22));
        MatrixXf M_6 = strassen((A_21-A_11), (B_11+B_12));
        MatrixXf M_7 = strassen((A_12-A_22), (B_21+B_22));

        C.block(0, 0, block_size, block_size) = M_1 + M_4 - M_5 + M_7;
        C.block(0, block_size, block_size, block_size) = M_3 + M_5;
        C.block(block_size, 0, block_size, block_size) = M_2 + M_4;
        C.block(block_size, block_size, block_size, block_size) = M_1 - M_2 + M_3 + M_6;
   }

	return C;
}

/*
 * f)
 *
 * O(N^log2(7))
 *
 */

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

