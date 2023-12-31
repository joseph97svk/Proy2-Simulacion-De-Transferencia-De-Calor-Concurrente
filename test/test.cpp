/* Copyright 2022 Joseph Stuart Valverde Kong.
Universidad de Costa Rica. CC BY 4.0 */
#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>

/**
 * @brief define a Matrix as a vector of a vector
 * 
 * @tparam dataType 
 */
template <typename dataType>
using Matrix = std::vector<std::vector<dataType>>;

/**
 * @brief returns address of matrix with
 * data extracted from file with file name
 * 
 * @param fileName of file where matrix is to be extracted
 * @return Matrix<double>* 
 */
Matrix<double>* getMatrix(std::string fileName);

/**
 * @brief compares two matrixes
 * @details returns false and reports on first error
 * 
 * @param matrixA first matrix
 * @param matrixB second matrix
 * @param epsilom acceptable error betwen both matrixes
 * @return true if both matrixes are equal withing parameters
 * @return false if an error was found
 */
bool compareMatrixes(Matrix<double>& matrixA,
    Matrix<double>& matrixB, double epsilom);

int main(int argc, char* argv[]) {
  // check if argument count is as should
  if (argc != 4) {
    std::cout << "Wrong amount of arguments given!" << std::endl;
    return EXIT_FAILURE;
  }

  try {
    // get file names
    std::string fileNameA = argv[1];
    std::string fileNameB = argv[2];

    // get epsilon
    double epsilon = std::stod(argv[3]);

    // extract matrixes from files
    Matrix<double> matrixA; 
    getMatrix(matrixA, fileNameA);
    Matrix<double> matrixB; 
    getMatrix(matrixB, fileNameB);

    // compare matrixes
    if (compareMatrixes(matrixA, matrixB, epsilon)) {
      std::cout << "All good!" << std::endl;
    }
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << std::endl;
  }
}

// returns address of matrix with data extracted from file with file name
void getMatrix(Matrix<double>& dataMatrix, std::string fileName) {
  // open binary file
  std::ifstream file;
  file.open(fileName, std::ios::binary);

  // check if binary file is opened
  if (!file.is_open()) {
    throw std::runtime_error("Binary file could not be opened for a given job");
  }

  size_t rowAmount = 0, colAmount = 0;

  // read matrix dimensions
  file.read(reinterpret_cast<char*>(&rowAmount), sizeof(size_t));
  file.read(reinterpret_cast<char*>(&colAmount), sizeof(size_t));

  // set size of matrix
  dataMatrix.resize(rowAmount);

  // for each row
  for (size_t row = 0; row < rowAmount; ++row) {
    // resize each row
    dataMatrix[row].resize(colAmount);

    // read from file and into matrix
    file.read(reinterpret_cast<char*>
        (&(dataMatrix[row][0])), sizeof(double) * colAmount);
  }

  // close file
  file.close();
}

// compares two matrixes and returns if they are equal within given parameters
bool compareMatrixes(Matrix<double>& matrixA,
    Matrix<double>& matrixB, double epsilom) {
  // check if row sizes coincide
  if (matrixA.size() != matrixB.size()) {
    std::cout << "Error: matrixes of mismatching row amount\n"
    << "matrix A size: " << matrixA.size() << std::endl
    << "matrix B size: " << matrixB.size() << std::endl;
    return false;
  }

  // check if col sizes coincide
  if (matrixA[0].size() != matrixB[0].size()) {
    std::cout << "Error: matrixes of mismatching row size\n"
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

      // if not within epsilon
      if (difference > epsilom) {
        // report failure
        std::cout << "Error: cell with wrong results at: "
        << row << ", " << col << std::endl;

        //  stop checking
        return false;
      }
    }
  }

  return true;
}
