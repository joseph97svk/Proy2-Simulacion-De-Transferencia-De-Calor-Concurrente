/* Copyright 2022 Joseph Stuart Valverde Kong.
Universidad de Costa Rica. CC BY 4.0 */
#include <cstdlib>
#include <utility>
#include <ctime>
#include <iomanip>

#include "TermalTransferSimulation.hpp"
#include "TermalData.hpp"

/**
 * @brief runs a stage of the simulation
 * 
 * @param data to be processed
 * @param newData where new data will be stored
 * @param jobInformation parameters for processing
 */
void runStage(Matrix<double>& data, Matrix<double>& newData,
JobInformation* jobInformation);

/**
 * @brief checks if equilibrium has been reached according to given parameters
 * 
 * @param data original data
 * @param newData new data
 * @param equilibriumPointSentivity epsilon of max change to be within equilibirum
 * @return true if in equilibrium
 * @return false if not within equilibrium
 */
bool checkEquilibrium(Matrix<double>& data, Matrix<double>& newData,
const double equilibriumPointSentivity);

/**
 * @brief writes report of all simulations
 * 
 * @param jobsInformation information on all jobs
 * @param fileName name of file where all jobs were listed
 */
void writeReport(std::vector<JobInformation>& jobsInformation, std::string fileName);

/**
 * @brief gets the time formated as string
 * 
 * @param seconds amount of seconds taken
 * @return std::string
 */
static std::string format_time(const time_t seconds);

// returns a vector of jobs found within file from file name
std::vector<JobInformation>* TermalTransferS::getJobData
(std::string& fileName) {
  // create vector for jobs
  std::vector<JobInformation>* dataVector = new std::vector<JobInformation>();

  // open file
  std::ifstream file;
  file.open(fileName);

  std::string extension;

  bool extensionNotFound = true;
  int position = fileName.size() - 1;

  // find where the extension ends
  while (extensionNotFound && position > 0) {
    if (fileName[position] == 47) {
      extensionNotFound = false;
    }

    position--;
  }

  // get the extension out of the file name
  if (position != 0) {
    position+=2;
    extension = fileName.substr(0, position);
  }

  int dataPosition = 0;
  int jobsFoundAmount  = 0;
  std::string currentData;

  // while data can be extracted from the file
  while (file >> currentData) {
    // place a JobInformation on the vector
    dataVector->emplace_back(JobInformation());

    // place data on the jobInformation
    switch (dataPosition) {
      // place the file name
      case 0:
        (*dataVector)[jobsFoundAmount].fileName = currentData.c_str();
        break;
      // place the time
      case 1:
        (*dataVector)[jobsFoundAmount].stageTimeDuration = stod(currentData);
        break;
      // place the termal difusivity
      case 2:
        (*dataVector)[jobsFoundAmount].termalDifusivity = stod(currentData);
        break;
      // place the cell dimensions
      case 3:
        (*dataVector)[jobsFoundAmount].cellDimensions = stod(currentData);
        break;
      // place epsilon
      case 4:
        (*dataVector)[jobsFoundAmount].
            equilibriumPointSentivity = stod(currentData);
        break;
      default:
        break;
    }

    dataPosition++;

    // advance to next job if all current data already placed
    if (dataPosition == 5) {
      jobsFoundAmount++;
      dataPosition = 0;

      // place the same extension on all JobInformations
      (*dataVector)[jobsFoundAmount].extension = extension;
    }
  }

  // shrink to size
  dataVector->resize(jobsFoundAmount);

  return dataVector;
}

// processess all jobs in the vector
void TermalTransferS::processAllJobs(std::vector<JobInformation>* jobs, std::string& fileName) {
  int jobsAmount = jobs->size();

  // for all jobs
  for (int currentJob = 0; currentJob < jobsAmount; ++currentJob) {
    // process them
    processJob(&(*jobs)[currentJob]);
  }

  writeReport(*jobs, fileName);
}

/**
 * @brief prints matrix on console
 * 
 * @param matrix to be printed
 */
void debugMatrix(Matrix<double>& matrix) {
  size_t rowAmount = matrix.size(),
      colAmount = matrix[0].size();

  std::cout << rowAmount << "," << colAmount << std::endl;

  for (size_t row = 0; row < rowAmount; ++row) {
    for (size_t col = 0; col < colAmount; ++col) {
      std::cout << matrix[row][col] << " ";
    }
    std::cout << std::endl;
  }

  std::cout << std::endl;
}

// processess the given job
void TermalTransferS::processJob(JobInformation* jobInformation) {
  // get the matrix out of the file
  Matrix<double>* data = getMatrix(jobInformation);

  // set a new matrix for new data
  Matrix<double>* newData = new Matrix<double>;

  newData->resize(data->size());

  // copy all data from previous to new (to account for edges)
  for (size_t row = 0; row < newData->size(); ++row) {
    (*newData)[row].resize((*data)[0].size());
    for (size_t col = 0; col < (*newData)[row].size(); ++col) {
      (*newData)[row][col] =
      (*data)[row][col];
    }
  }

  bool inEquilibrium = false;

  int stageCount = 0;

  // while the matrix is not in equilibrium
  while (!inEquilibrium) {
    // run a stage/find state
    runStage(*data, *newData, jobInformation);

    // check if in equilibrium
    inEquilibrium =
    checkEquilibrium(*data, *newData,
    jobInformation->equilibriumPointSentivity);

    // if not in equilibrium
    if (!inEquilibrium) {
      /* swap matrixes (data will always have old data)
       * old data will be overwritten and lost within the matrixes
      */
      std::swap(data, newData);
    }

    // increase stages processed
    stageCount++;
  }

  jobInformation->stateAmountRequired = stageCount;

  // write results on new file
  writeMatrixOnFile(*newData, jobInformation);
  delete(data);
}

// returns a matrix of the data found in the file data
void runStage(Matrix<double>& data, Matrix<double>& newData,
JobInformation* jobInformation) {
  // get sizes
  size_t rowAmount = data.size(),
      colAmount = data[0].size();

  // get formula parameters
  double time = jobInformation->stageTimeDuration,
      alpha = jobInformation->termalDifusivity,
      dimensions =
      jobInformation->cellDimensions * jobInformation->cellDimensions;

  // for each row
  for (size_t row = 1; row < rowAmount - 1; ++row) {
    // TODO(me): run a thread for each of these
    // for each column
    for (size_t col = 1; col < colAmount - 1; ++col) {
      // find the sum of neighbors
      double sumOfNeightbors =
      data[row][col - 1] + data[row - 1][col] +
      data[row][col + 1] + data[row + 1][col] -
      (4 * (data[row][col]));

      // find the new state of the block
      newData[row][col] = data[row][col] +
      (((time * alpha)/dimensions) * sumOfNeightbors);
    }
  }
}

// writes matrix data on binary file
bool checkEquilibrium(Matrix<double>& data, Matrix<double>& newData,
const double equilibriumPointSentivity) {
  // get sizes
  size_t rowAmount = data.size(),
      colAmount = data[0].size();

  // for each row
  for (size_t row = 1; row < rowAmount - 1; ++row) {
    // for each column
    for (size_t col = 1; col < colAmount - 1; ++col) {
      // find the difference
      double difference = data[row][col] - newData[row][col];
      // ensure it is positive
      if (difference < 0) {
        difference = -difference;
      }
      // check if it is within epsilom
      if (difference > equilibriumPointSentivity) {
        return false;
      }
    }
  }
  // always within epsilom return true
  return true;
}

// returns a matrix of the data found in the file data
Matrix<double>* TermalTransferS::getMatrix(JobInformation* jobInformation) {
  // get complete path with file name
  std::string fileName = "jobs/" + jobInformation->fileName;

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

// writes matrix data on binary file
void TermalTransferS::writeMatrixOnFile(Matrix<double>& dataMatrix,
    JobInformation* jobInformation) {
  // get dimensions
  size_t rowAmount = dataMatrix.size(),
      colAmount = dataMatrix[0].size();

  // get new file name
  std::string newFileName = jobInformation->extension +
      jobInformation->fileName.substr(0, jobInformation->fileName.size() - 4)
      + "-" + std::to_string(jobInformation->stateAmountRequired) +".bin";

  // create binary file
  std::ofstream newFile(newFileName, std::ios::binary);

  // write matrix dimensions
  newFile.write(reinterpret_cast<const char*>(&rowAmount), sizeof(size_t));
  newFile.write(reinterpret_cast<const char*>(&colAmount), sizeof(size_t));

  // ensure size is correct
  dataMatrix.resize(rowAmount);

  // for each row
  for (size_t row = 0; row < rowAmount; ++row) {
    // write all columns
    newFile.write(reinterpret_cast<const char*>
        (&dataMatrix[row][0]), (sizeof(double) * (colAmount)));
  }

  // delete matrix
  delete(&dataMatrix);

  // close file
  newFile.close();
}

// erases all data related to JobInformation struct
void TermalTransferS::eraseJobData(std::vector<JobInformation>* jobData) {
  // erases job data
  delete(jobData);
}

void writeReport(std::vector<JobInformation>& jobsInformation, std::string fileName) {
  // open file
  std::ifstream file;
  file.open(fileName);

  std::string newFileName = jobsInformation[0].extension +
      fileName.substr(0, fileName.size() - 4) +".tsv";;

  // create binary file
  std::ofstream newFile(newFileName, std::ios::binary);

  std::string current;

  int job = 0;
  // for each line read
  while (getline(file, current)) {
    // write it to the new file
    newFile << current;
    // get the time
    double totalTime = jobsInformation[job].stageTimeDuration *
    jobsInformation[job].stateAmountRequired;
    // add time information
    newFile << std::setfill(' ') << std::setw(30) << format_time(totalTime);
    // next line
    newFile << "\n";
    ++job;
  }
  file.close();
  newFile.close();
}

static std::string format_time(const time_t seconds) {
  // TODO(any): Using C until C++20 std::format() is implemented by compilers
  char text[48];  // YYYY/MM/DD hh:mm:ss
  const std::tm& gmt = * std::gmtime(&seconds);
  snprintf(text, sizeof(text), "%04d/%02d/%02d\t%02d:%02d:%02d", gmt.tm_year
    - 70, gmt.tm_mon, gmt.tm_mday - 1, gmt.tm_hour, gmt.tm_min, gmt.tm_sec);
  return text;
}
