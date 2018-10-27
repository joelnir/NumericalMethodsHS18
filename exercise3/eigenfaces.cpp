#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <eigen3/Eigen/Dense>

using namespace Eigen;

VectorXd load_pgm(const std::string &filename) {
	// returns a picture as a flattened vector

	int row = 0, col = 0, rows = 0, cols = 0;

	std::ifstream infile(filename);
	std::stringstream ss;
	std::string inputLine = "";

	// First line : version
	std::getline(infile,inputLine);

	// Second line : comment
	std::getline(infile,inputLine);

	// Continue with a stringstream
	ss << infile.rdbuf();
	// Third line : size
	ss >> cols >> rows;

	VectorXd picture(rows*cols);

	// Following lines : data
	for(row = 0; row < rows; ++row) {
		for (col = 0; col < cols; ++col) {
			int val;
			ss >> val;
			picture(col*rows + row) = val;
		}
	}

	infile.close();
	return picture;
}

int main() {

	int h = 231;
	int w = 195;
	int M = 15;

	MatrixXd faces(h*w, M);
	MatrixXd A(h*w, M);

	VectorXd meanFace(h*w);
	VectorXd sumFace = VectorXd::Zero(h*w);

	// loads pictures as flattened vectors into faces
	for (int i=0; i<M; i++) {
		std::string filename = "./basePictures/subject"+
			std::to_string(i+1) + ".pgm";
		VectorXd flatPic = load_pgm(filename);
		faces.col(i) = flatPic;

        sumFace = sumFace += flatPic;
	}

    // Compute mean
    meanFace = (1/M) * sumFace;

    for(int i = 0; i < M; ++i){
        A.col(i) = faces.col(i) - meanFace;
    }

    JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);
    MatrixXd eigenFaces = svd.matrixU();

	// try to recognize a test face
    std::string testPicName = "./testPictures/Narutowicz.pgm";
	VectorXd newFace = load_pgm(testPicName);

    VectorXd proj(M);
    proj = eigenFaces.transpose() * (newFace - meanFace);

    int minK = -1;
    double minDist = 0;
    double dist;
    VectorXd proj_k(M);;
    for(int k = 0; k < M; ++k){
        proj_k =  eigenFaces.transpose() * A.col(k);

        dist = (proj-proj_k).norm();

        if((minK == -1) || (dist < minDist)){
            minK = k+1;
            minDist = dist;
        }
    }

    std::cout << "Minimal k: " << minK << std::endl;
}
