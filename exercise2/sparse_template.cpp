//#include <chrono>
#include <functional>
#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

std::vector<Eigen::Triplet<double>> MakeTripletList(int n) {
	int nnz = 3 * n - 2;
	std::vector<Eigen::Triplet<double>> tripletList(nnz);

	// TODO: Task (b)
	// ...
	// ...
	// ...
	// ...
	// ...

	return tripletList;
}

double Runtime(const std::function<void(void)> &f) {

	// TODO: Task (c)
	// ...
	// ...
	// ...
	// ...
	// ...
	// ...
	// ...
	// ...
	// ...
	// ...

    return .0;	// dummy return value
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
		// TODO: Task (c)
		// ...
		// ...
		// ...

		// dense LU decomposition
		Eigen::MatrixXd D(S);
		// TODO: Task (c)
		// ...
		// ...
		// ...

		// benchmark
		// TODO: Task (c)
		// ...
	}

	std::cout << "Runtime in seconds using storage format..." << std::endl;
	std::cout << "...sparse: " << runtimeSparse << std::endl;
	std::cout << "...dense:  " << runtimeDense << std::endl;

	return 0;
}
