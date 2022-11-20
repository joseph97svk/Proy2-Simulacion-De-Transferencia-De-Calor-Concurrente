#include <cstdlib>

#include "TermalTransferSimulation.hpp"
#include "TermalData.hpp"

void runStage (Matrix<double>& data, Matrix<double>& newData,
JobInformation* jobInformation);

bool checkEquilibrium(Matrix<double>& data, Matrix<double>& newData,
const double equilibriumPointSentivity);

std::vector<JobInformation*>* TermalTransferS::getJobData(std::string& fileName) {
  std::vector<JobInformation*>* dataVector = new std::vector<JobInformation*>();
  
  std::ifstream file;
  file.open(fileName);

  std::string currentData;
  int dataPosition = 0;
  int jobsFoundAmount  = 0;

  std::string extension;

  bool extensionNotFound = true;
  int position = fileName.size() - 1;
  while (extensionNotFound && position > 0) {
    if (fileName[position] == 47) {
      extensionNotFound = false;
    }

    position--;
  }

  if (position != 0) {
    position+=2;
    extension = fileName.substr(0, position);
  }

  while (file >> currentData) {
    dataVector->push_back(new JobInformation());
    
    switch (dataPosition) {
      case 0:
        (*dataVector)[jobsFoundAmount]->fileName = currentData.c_str();
        break;
      case 1:
        (*dataVector)[jobsFoundAmount]->stageTimeDuration = stod(currentData);
        break;
      case 2:
        (*dataVector)[jobsFoundAmount]->termalDifusivity = stod(currentData);
        break;
      case 3:
        (*dataVector)[jobsFoundAmount]->cellDimensions = stod(currentData);
        break;
      case 4:
        (*dataVector)[jobsFoundAmount]->equilibriumPointSentivity = stod(currentData);
        break;
      default:
        break;
    }

    (*dataVector)[jobsFoundAmount]->extension = extension;
    
    dataPosition++;
    if (dataPosition == 5) {
      jobsFoundAmount++;
      dataPosition = 0;
    }
  }

  dataVector->resize(jobsFoundAmount);
  
  return dataVector;
}

void TermalTransferS::processAllJobs(std::vector<JobInformation*>* jobs) {
  int jobsAmount = jobs->size();

  std::cout << jobsAmount << std::endl;

  for (int currentJob = 0; currentJob < jobsAmount; ++currentJob) {
    processJob((*jobs)[currentJob]);
  }
}

void TermalTransferS::processJob(JobInformation* jobInformation) {
  Matrix<double>* data = getMatrix(jobInformation);
  
  Matrix<double>* newData = new Matrix<double>;

  newData->resize(data->size());
  for (size_t row = 0; row < newData->size(); ++row) {
    (*newData)[row].resize((*data)[0].size());
  }

  bool inEquilibrium = false;

  int stageCount = 0;

  while (!inEquilibrium) {
    runStage(*data, *newData, jobInformation);

    inEquilibrium =
    checkEquilibrium(*data, *newData,
    jobInformation->equilibriumPointSentivity);

    Matrix<double>* temp = newData;
    newData = data;
    data = temp;
    stageCount++;
  }

  if (stageCount % 2 == 0) {
    Matrix<double>* temp = newData;
    newData = data;
    data = temp;
  }
  std::cout << "job done" << std::endl;
  writeMatrixOnFile(*newData, jobInformation);
}

void runStage (Matrix<double>& data, Matrix<double>& newData,
JobInformation* jobInformation) {
  size_t rowAmount = data.size(),
      colAmount = data[0].size();

  double time = jobInformation->stageTimeDuration,
      alpha = jobInformation->termalDifusivity,
      dimensions =
      jobInformation->cellDimensions * jobInformation->cellDimensions;


  for (size_t row = 1; row < rowAmount - 1; ++row) {
    // run a thread for each of these
    for (size_t col = 1; col < colAmount - 1; ++col) {
      double sumOfNeightbors =
      data[row][col - 1] + data[row - 1][col] + data[row][col + 1] + data[row + 1][col];

      newData[row][col] = data[row][col] + 
      (((time * alpha)/dimensions) * sumOfNeightbors);
    }
  }
}

bool checkEquilibrium(Matrix<double>& data, Matrix<double>& newData,
const double equilibriumPointSentivity) {
  size_t rowAmount = data.size(),
      colAmount = data[0].size();

  for (size_t row = 0; row < rowAmount; ++row) {
    for (size_t col = 0; col < colAmount; ++col) {
      if (abs(data[row][col] - newData[row][col])
      > equilibriumPointSentivity) {
        return false;
      }
    }
  }

  return true;
}

Matrix<double>* TermalTransferS::getMatrix(JobInformation* jobInformation) {
  std::string fileName = "jobs/" + jobInformation->fileName;

  std::ifstream file;
  file.open(fileName, std::ios::binary);

  size_t rowAmount = 0, colAmount = 0;

  file.read(reinterpret_cast<char*>(&rowAmount), sizeof(size_t));
  file.read(reinterpret_cast<char*>(&colAmount), sizeof(size_t));

  Matrix<double>* dataMatrix = new Matrix<double>;

  Matrix<double>& matrixRef = *dataMatrix;

  std::cout << ">" << rowAmount << std::endl;

  matrixRef.resize(rowAmount);

  std::cout << ">" << matrixRef.size() << std::endl;
  
  for (size_t row = 0; row < rowAmount; ++row) {
    matrixRef[row].resize(colAmount);
    for (size_t col = 0; col < colAmount; ++col) {
      file.read(reinterpret_cast<char*>(&(matrixRef[row][col])), sizeof(double));
    }
  } 

  return dataMatrix;
}

void TermalTransferS::writeMatrixOnFile(Matrix<double>& dataMatrix, JobInformation* jobInformation) {
  size_t rowAmount = dataMatrix.size(),
      colAmount = dataMatrix[0].size();

  std::string newFileName =
      jobInformation->fileName.substr(0, jobInformation->fileName.size() - 4) + "-k.bin";

  std::ofstream newFile(newFileName, std::ios::binary);

  newFile.write(reinterpret_cast<char*>(&rowAmount), sizeof(size_t));
  newFile.write(reinterpret_cast<char*>(&colAmount), sizeof(size_t));

  for (size_t row = 0; row < rowAmount; ++row) {
    newFile.write(reinterpret_cast<char*>(&dataMatrix[row]), sizeof(double) * colAmount);
  }
}