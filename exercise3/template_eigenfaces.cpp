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
	VectorXd meanFace(h*w);
	
	// loads pictures as flattened vectors into faces
	for (int i=0; i<M; i++) {
		std::string filename = "./basePictures/subject"+ 
			std::to_string(i+1) + ".pgm";
		VectorXd flatPic = load_pgm(filename);
		faces.col(i) = flatPic;
		
		// TODO: Point (b)
	}

	
	// TODO: Point (e)

	// try to recognize a test face
	string testPicName = "./testPictures/subject01.happy.pgm";
	VectorXd newFace = load_pgm(testPicName);

	// TODO: Point (f)
	
	// TODO: Point (g)

}
