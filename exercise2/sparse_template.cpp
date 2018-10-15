#include <chrono>
#include <functional>
#include <iostream>
#include <vector>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>

/*
 * a)
 *
 * (Assume n >= 2)
 * (3*(n - 2) + 4)*2 + (n + 1) = 7n - 3
 */

std::vector<Eigen::Triplet<double>> MakeTripletList(int n) {
	int nnz = 3 * n - 2;
	std::vector<Eigen::Triplet<double>> tripletList(nnz);

    for(int i = 0; i < n; ++i){
        if(i > 0){
            tripletList.push_back(Eigen::Triplet<double>(i, i-1, -1.0));
        }

        tripletList.push_back(Eigen::Triplet<double>(i, i, 2.0));

        if(i < (n-1)){
            tripletList.push_back(Eigen::Triplet<double>(i, i+1, -1.0));
        }
    }

	return tripletList;
}

double Runtime(const std::function<void(void)> &f) {
    auto start = std::chrono::system_clock::now();

    // Run function
    f();

    auto end = std::chrono::system_clock::now();
    auto elapsed = end - start;
    return elapsed.count();	// dummy return value
}

template <class T>
std::ostream & operator<< (std::ostream &os, const std::vector<T> &v) {
	os << "[";
	if (!v.empty()) {
		os << v[0];
		for (int i = 1; i < v.size(); ++i) os << ", " << v[i];
	}
    os << "]";

    return os;
}

int main() {
	// print small example of the tridiagonal matrix
	int m = 4;
	std::vector<Eigen::Triplet<double>> tripletList = MakeTripletList(m);
	Eigen::SparseMatrix<double> S_(m, m);
	S_.setFromTriplets(tripletList.begin(), tripletList.end());
	std::cout << "If n = " << m << ", then T equals" << std::endl;
	std::cout << Eigen::MatrixXd(S_) << std::endl;

	// matrix sizes for benchmark
	std::vector<int> N = {64, 128, 256, 512};
	std::cout << "LU decomposition of T, where n = " << N << std::endl;

	// set up variables for runtime measurement
	std::vector<double> runtimeSparse;
	std::vector<double> runtimeDense;

	for (int n : N) {
		tripletList = MakeTripletList(n);

		// sparse LU decomposition
		Eigen::SparseMatrix<double> S(n, n);
		S.setFromTriplets(tripletList.begin(), tripletList.end());

		// dense LU decomposition
		Eigen::MatrixXd D(S);

        Eigen::SparseLU<Eigen::SparseMatrix<double> > sparseLU;
        Eigen::FullPivLU<Eigen::MatrixXd> denseLU;

		// benchmark
        auto sparseFunc = [&S, &sparseLU](){
            sparseLU.compute(S);
        };
        auto denseFunc = [&D, &denseLU](){
            denseLU.compute(D);
        };

        runtimeSparse.push_back(Runtime(sparseFunc));
        runtimeDense.push_back(Runtime(denseFunc));
	}

	std::cout << "Runtime in seconds using storage format..." << std::endl;
	std::cout << "...sparse: " << runtimeSparse << std::endl;
	std::cout << "...dense:  " << runtimeDense << std::endl;

	return 0;
}

/*
 * d)
 *
 * Considering nnz in= O(n) =>
 *
 * Dense LU: O(n^3)
 * Sparse LU: ~O(nnz^2) = O(n^2)
 *
 */
