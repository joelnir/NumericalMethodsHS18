#include <iostream>

#include <Eigen/Dense>

/* @brief 
 * @param[in] z An $n$-dimensional vector containing one side of input data
 * @param[in] c An $n$-dimensional vector containing the other side of input data
 * @param[out] x The vector of parameters $(\alpha,\beta)$, intercept and slope of the line fitted
 */
Eigen::Vector2d lsqEst(const Eigen::VectorXd &z, const Eigen::VectorXd &c) {
	Eigen::Vector2d x;

	// TODO: Task (b)
	// ...
	// ...
	// ...

	return x;
}

int main() {
    int n = 10;
    Eigen::VectorXd z(n), c(n);
    for(int i = 0; i < n; ++i) {
		z(i) = i + 1;
		c(i) = n - i;
	}

	Eigen::Vector2d x = lsqEst(z, c);

	std::cout << "alpha = " << x(0) << std::endl;
	std::cout << "beta = "  << x(1) << std::endl;

	return 0;
}
