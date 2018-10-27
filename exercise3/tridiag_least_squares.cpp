#include <iostream>

#include <eigen3/Eigen/Dense>

/* @brief
 * @param[in] z An $n$-dimensional vector containing one side of input data
 * @param[in] c An $n$-dimensional vector containing the other side of input data
 * @param[out] x The vector of parameters $(\alpha,\beta)$, intercept and slope of the line fitted
 */
Eigen::Vector2d lsqEst(const Eigen::VectorXd &z, const Eigen::VectorXd &c) {
	Eigen::Vector2d x;
    int n = z.size();
    Eigen::MatrixXd A(n, 2);

    Eigen::VectorXd z_hat1(n);
    Eigen::VectorXd z_hat2(n);

    z_hat1 << 0, z.head(n-1);
    z_hat2 << z.tail(n-1), 0;

    A << z, (z_hat1 + z_hat2);

    // Solve LSQ Ax - c
    x = (A.transpose() * A).llt().solve(A.transpose() * c);
	return x;
}

int main() {
    int n = 3;
    Eigen::VectorXd z(n), c(n);
    for(int i = 0; i < n; ++i) {
		z(i) = i + 1;
		c(i) = n - i;
	}

    std::cout << "z:" << z << std::endl;
    std::cout << "c:" << c << std::endl;

	Eigen::Vector2d x = lsqEst(z, c);

	std::cout << "alpha = " << x(0) << std::endl;
	std::cout << "beta = "  << x(1) << std::endl;

	return 0;
}
