#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>

template <typename dataType>
using Matrix = std::vector<std::vector<dataType>>;

Matrix<double>* getMatrix(std::string fileName);

bool compareMatrixes(Matrix<double>& matrixA,
    Matrix<double>& matrixB, double epsilom);

int main (int argc, char* argv[]) {
  if (argc != 4) {
    std::cout << "Wrong amount of arguments given!" << std::endl;
    return EXIT_FAILURE;
  }

  // get file names
  std::string fileNameA = argv[1];
  std::string fileNameB = argv[2];

  double epsilon = std::stod(argv[3]);

  // extract matrixes from files
  Matrix<double>* matrixA = getMatrix(fileNameA);
  Matrix<double>* matrixB = getMatrix(fileNameB);
  
  // compare matrixes
  if (compareMatrixes(*matrixA, *matrixB, epsilon)) {
    std::cout << "All good!" << std::endl;
  }

  delete matrixA;
  delete matrixB;
}

Matrix<double>* getMatrix(std::string fileName) {
  // open binary file
  std::ifstream file;
  file.open(fileName, std::ios::binary);

  size_t rowAmount = 0, colAmount = 0;

  // read matrix dimensions
  file.read(reinterpret_cast<char*>(&rowAmount), sizeof(size_t));
  file.read(reinterpret_cast<char*>(&colAmount), sizeof(size_t));

  // make new matrix
  Matrix<double>* dataMatrix = new Matrix<double>;

  // set reference for ease of use
  Matrix<double>& matrixRef = *dataMatrix;

  // set size of matrix
  dataMatrix->resize(rowAmount);

  // for each row
  for (size_t row = 0; row < rowAmount; ++row) {
    // resize each row
    matrixRef[row].resize(colAmount);
    // for each column
    for (size_t col = 0; col < colAmount; ++col) {
      // read from file and into matrix
      file.read(reinterpret_cast<char*>
          (&(matrixRef[row][col])), sizeof(double));
    }
  }

  // close file
  file.close();

  return dataMatrix;
}

bool compareMatrixes(Matrix<double>& matrixA,
    Matrix<double>& matrixB, double epsilom) {
  if (matrixA.size() != matrixB.size()) {
    std::cout << "Error: matrixes of mismatching row amount"
    << "matrix A size: " << matrixA.size() << std::endl
    << "matrix B size: " << matrixB.size() << std::endl;
    return false;
  }

  if (matrixA[0].size() != matrixB[0].size()) {
    std::cout << "Error: matrixes of mismatching row size"
    << "matrix A size: " << matrixA[0].size() << std::endl
    << "matrix B size: " << matrixB[0].size() << std::endl;
    return false;
  }
  
  size_t rowAmount = matrixA.size(),
    colAmount = matrixA[0].size();

  for (size_t row = 0; row < rowAmount; ++row) {
    // for each column
    for (size_t col = 0; col < colAmount; ++col) {
      double difference = matrixA[row][col] - matrixB[row][col];
      // ensure it is positive
      if (difference < 0) {
        difference = -difference;
      }

      if (difference > epsilom) {
        std::cout << "Error: cell with wrong results at: "
        << row << ", " << col << std::endl;

        return false;
      }
    }
  }

  return true;
}