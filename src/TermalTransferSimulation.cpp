/* Copyright 2022 Joseph Stuart Valverde Kong.
Universidad de Costa Rica. CC BY 4.0 */
#include <cstdlib>
#include <utility>
#include <ctime>
#include <iomanip>
#include <omp.h>  // NOLINT
#include <thread>
#include <mpi.h>  // NOLINT
#include <stdexcept>

#include "TermalTransferSimulation.hpp"
#include "TermalData.hpp"

#define CAST(type, unit)\
  reinterpret_cast<type>(unit)

/**
 * @brief processess all jobs in the vector with one thread
 * 
 * @param jobs vector of jobs to be processed
 * @param fileName name of file where all jobs/plates are located
 * @param threadAmount amount of threads expected to run
 */
void processAllJobsST(std::vector<JobInformation>& jobs,
std::string& fileName, int32_t threadAmount);

/**
 * @brief processess the given job
 * 
 * @param jobInformation job to be processed
 * @param threadAmount amount of threads expected to run
 */
void processJobST(JobInformation& jobInformation, int32_t threadAmount);

/**
 * @brief administrative process handling task issue and process coordination
 * 
 * @param fileName name of file with all job/plate information
 * @param size amount of processes in program
 */
void mpiRank0(std::string fileName, const int32_t size);

/**
 * @brief sends data to all processes that request it
 * 
 * @param jobData vector of all jobs
 * @param size amount of processes running for the program
 */
void sendData(std::vector<JobInformation>& jobData, const int32_t size);

/**
 * @brief sends stop conditions to stop all processes
 * @details stop condition is time -1
 * @param size amount of processes
 */
void stopProcessess(const int32_t size);

/**
 * @brief worker processes assigned to processing
 * 
 * @param rank process number or ID
 * @param threadAmount amount of desired threads
 * @param fileName name of file with all plates/jobs
 */
void mpiRankAny(const int32_t rank,
    const int32_t threadAmount, std::string fileName);

/**
 * @brief Receives data for worker process from administrating process
 * 
 * @param stages reference to amount of stages already processed
 * @return int32_t 
 */
int32_t receiveData(size_t& stages);

/**
 * @brief Gets the information for a job at a given position
 * 
 * @param job job where all information will be stored
 * @param fileName name of file with all jobs/plates
 * @param jobPosition position of job/plate within file
 */
void getSingleJobData(JobInformation& job,
    std::string& fileName, size_t jobPosition);

/**
 * @brief Gets the path out of a file name
 * 
 * @param fileName complete file name with path
 * @return std::string path of the file
 */
std::string getFilePath(std::string fileName);

/**
 * @brief runs a stage of the simulation and
 * checks how many cell changes break equilibrium
 * 
 * @param data to be processed
 * @param newData where new data will be stored
 * @param jobInformation parameters for processing
 * @param equilibriumPointSentivity epsilon of change that must not be exceeded
 * @return int amount of cells breaking equilibrium
 */
int runStage(Matrix<double>& data, Matrix<double>& newData,
    JobInformation& jobInformation, const double equilibriumPointSentivity);

/**
 * @brief writes report of all simulations
 * 
 * @param jobsInformation information on all jobs
 * @param fileName name of file where all jobs were listed
 */
void writeReport(std::vector<JobInformation>& jobsInformation,
    std::string fileName);

/**
 * @brief 
 * 
 * @tparam datatype 
 * @param dataA 
 * @param dataB 
 */
template<typename datatype>
void mySwap(datatype& dataA, datatype& dataB) {
  datatype& temp = dataA;
  dataA = dataB;
  dataB = temp;
}

// runs termal transfer simulation
void TermalTransferS::runTermalTransferSimulation(int& argc, char**& argv) {
  if (argc > 1) {
    // get file name
    std::string fileName(argv[1]);

    // set default amount of threads
    int32_t threadAmount = std::thread::hardware_concurrency();

    // if thread is given
    if (argc == 3) {
      // then change it
      threadAmount = std::stoi(argv[2]);
    }

    // Initiate mpi connection environment
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
      throw std::runtime_error("Could not initialize MPI enviroment");
    }

    int rank = 0, size = 0;

    // get the rank of the current process
    if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != EXIT_SUCCESS) {
      throw std::runtime_error("Could not get process rank");
    }

    // get the amount of processes
    if (MPI_Comm_size(MPI_COMM_WORLD, &size) != EXIT_SUCCESS) {
      throw std::runtime_error("Could not get process amount");
    }

    // if there is just one process
    if (size == 1) {
      std::vector<JobInformation> jobData;

      // get job data
      TermalTransferS::getJobData(jobData, fileName);

      // process everything
      TermalTransferS::processAllJobs(jobData, fileName, threadAmount);

    // if there are multiple process
    } else {
      // if rank 0
      if (rank == 0) {
        // (administrating process)
        mpiRank0(fileName, size);
      // for all other process
      } else {
        // (working process)
        mpiRankAny(rank, threadAmount, fileName);
      }
    }

    // finalize all process
    MPI_Finalize();
  } else {
    throw std::runtime_error
        ("No txt file with plate to process information given!");
  }
}

// administrative process handling task issue and process coordination
void mpiRank0(std::string fileName, const int32_t size) {
  // get jobs
  std::vector<JobInformation> jobData;
  TermalTransferS::getJobData(jobData, fileName);

  // send data to all process
  sendData(jobData, size);

  // stop processess with stop conditions
  stopProcessess(size);

  // write report
  writeReport(jobData, fileName);
}

// sends data to all processes that request it
void sendData(std::vector<JobInformation>& jobData, const int32_t size) {
  size_t time = 0;
  int32_t rankToSend = -1;

  // create array holding which job is being processed by which process
  int64_t kJobLocationArr[size]; // NOLINT

  int64_t jobAmount = (int64_t) jobData.size(); // NOLINT

  // for all jobs
  for (int64_t currentJob = 0; currentJob < jobAmount; ++currentJob) {
    // receive signal from processes ready to process
    MPI_Status status;
    MPI_Recv(&time, 1,
    MPI_LONG_LONG, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

    // save which process did which job
    kJobLocationArr[rankToSend] = currentJob;

    // set received time in correct job
    if (time != INT64_MAX) {
      jobData[kJobLocationArr[rankToSend] - 1].stateAmountRequired = time;
    }

    rankToSend = status.MPI_SOURCE;

    MPI_Send(&currentJob, /* amount */1, MPI_LONG_LONG,
        rankToSend, 1, MPI_COMM_WORLD);
  }
}

// sends stop conditions to stop all processes
void stopProcessess(const int32_t size) {
  int32_t rankToSend = -1;
  int64_t condition = -1;
  int64_t buffer = 0;

  // for all processess (stop conditions)
  for (int process = 0; process < size - 1; ++process) {
    MPI_Status status;

    // receive ready to process
    MPI_Recv(&buffer, 1,
    MPI_LONG_LONG, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

    rankToSend = status.MPI_SOURCE;

    // signal there is no more to do
    MPI_Send(&condition, /* amount */1, MPI_LONG_LONG,
        rankToSend, 1, MPI_COMM_WORLD);
  }
}

// worker processes assigned to processing
void mpiRankAny(const int32_t rank,
    const int32_t threadAmount, std::string fileName) {
  (void) rank;

  JobInformation job;
  int32_t processedJobs = 0;
  size_t stages = INT64_MAX;

  // while there are jobs
  while (true) {
    // receive all data for the job
    int64_t position =
    receiveData(stages);

    // end loop if no more data has been received
    if (position == -1) {
      break;
    }

    // get the job data for the given job
    getSingleJobData(job, fileName, position);

    // process job
    TermalTransferS::processJob(job, threadAmount);

    // get the amount of stages processed
    stages = job.stateAmountRequired;

    processedJobs++;
  }
}

// Receives data for worker process from administrating process
int32_t receiveData(size_t& stages) {
  int64_t jobPosition = -1;

  // send ready to process
  MPI_Send(&stages, 1, MPI_LONG_LONG, /*destination*/ 0, 0, MPI_COMM_WORLD);

  // receive position of the job/plate
  MPI_Recv(&jobPosition, 1,
    MPI_LONG_LONG, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  return jobPosition;
}

// returns a vector of jobs found within file from file name
void TermalTransferS::getJobData
    (std::vector<JobInformation>& dataVector, std::string& fileName) {
  // open file
  std::ifstream file;
  file.open(fileName);

  if (!file.is_open()) {
    throw std::runtime_error("Could not open txt file with job info");
  }

  // get file path
  std::string extension = getFilePath(fileName);

  int jobsFoundAmount  = 0;
  // read all jobs from the file
  do {
    // emplace back a new JobInformation struct
    dataVector.emplace_back(JobInformation());

    // add the extension (all files share the extension)
    dataVector[jobsFoundAmount].extension = extension;
  // continue while a job can be read
  } while (file >> dataVector[jobsFoundAmount++]);

  // adjust amount of jobs found
  jobsFoundAmount--;

  // shrink to size
  dataVector.resize(jobsFoundAmount);
}

// Gets the information for a job at a given position
void getSingleJobData(JobInformation& job,
    std::string& fileName, size_t jobPosition) {
  // open file
  std::ifstream file;
  file.open(fileName);

  // check if file is open
  if (!file.is_open()) {
    throw std::runtime_error("Could not open txt file with job info");
  }

  // get file path
  std::string extension = getFilePath(fileName);

  size_t position = 0;
  // read all jobs until finding one in position
  while (file >> job && position++ < jobPosition) {
  }

  // add the extension (all files share the extension)
  job.extension = extension;
}

// Gets the path out of a file name
std::string getFilePath(std::string fileName) {
  std::string extension;

  bool extensionNotFound = true;
  int position = fileName.size() - 1;

  // find where the extension ends
  while (extensionNotFound && position > 0) {
    if (fileName[position] == '/') {
      extensionNotFound = false;
    }

    position--;
  }

  // get the extension out of the file name
  if (position != 0) {
    position+=2;
    extension = fileName.substr(0, position);
  }

  return extension;
}

// processess all jobs in the vector
void TermalTransferS::processAllJobs(std::vector<JobInformation>& jobs,
    std::string& fileName, int32_t threadAmount) {
  int jobsAmount = jobs.size();

  if (threadAmount == 1) {
    // for all jobs
    for (int currentJob = 0; currentJob < jobsAmount; ++currentJob) {
      // process them
      processJobST(jobs[currentJob], threadAmount);
    }
    return;
  }

  // for all jobs
  for (int currentJob = 0; currentJob < jobsAmount; ++currentJob) {
    // process them
    processJob(jobs[currentJob], threadAmount);
  }

  // write tsv report
  writeReport(jobs, fileName);
}

// processess the given job
void TermalTransferS::processJob
    (JobInformation& jobInformation, int32_t threadAmount) {
  // get the matrix out of the file
  Matrix<double> data;
  getMatrix(data, jobInformation);

  // set a new matrix for new data
  Matrix<double> newData;
  newData.resize(data.size());

  size_t stageCount = 0, wrongAmount = 0;

  Matrix<double>* dataPointer = &data;
  Matrix<double>* newDataPointer = &newData;

  #pragma omp parallel num_threads(threadAmount) \
      default(none) shared(data, newData, dataPointer, newDataPointer, \
      jobInformation, stageCount, wrongAmount, std::cout)
  {  // NOLINT
    // copy all data from previous to new (to account for edges)
    #pragma omp for
    for (size_t row = 0; row < newData.size(); ++row) {
      newData[row].resize(data[0].size());

      // each thread takes a row at a time
      for (size_t col = 0; col < newData[row].size(); ++col) {
        newData[row][col] = data[row][col];
      }
    }

    bool localEquilibrium = false;

    // while the matrix is not in equilibrium
    while (!localEquilibrium) {
      // run a stage and receive wrong cell amount
      int localWrong = runStage(*dataPointer, *newDataPointer, jobInformation,
          jobInformation.equilibriumPointSentivity);

      // aggregate all non-equilibrium cells from all threads
      #pragma omp atomic
      wrongAmount += localWrong;

      // wait for all threads to report on their non-equilibrium cells
      #pragma omp barrier

      // check if in equilibrium
      localEquilibrium = wrongAmount == 0;

      // if not in equilibrium
      #pragma omp single
      if (!localEquilibrium) {
        /* swap matrixes (data will always have old data)
        * old data will be overwritten and lost within the matrixes
        */
        std::swap(*dataPointer, *newDataPointer);
      }

      // increase stages processed and reset coincidence amount
      #pragma omp single
      {
        stageCount++;
        wrongAmount = 0;
      }
    }
  }

  jobInformation.stateAmountRequired = stageCount;

  // write results on new file
  writeMatrixOnFile(newData, jobInformation);
}

// processess the given job
void processJobST(JobInformation& jobInformation, int32_t threadAmount) {
  (void) threadAmount;
  // get the matrix out of the file
  Matrix<double> data;
  TermalTransferS::getMatrix(data, jobInformation);

  // set a new matrix for new data
  Matrix<double> newData;
  newData.resize(data.size());

  size_t stageCount = 0, wrongAmount = 0;

  Matrix<double>* dataPointer = &data;
  Matrix<double>* newDataPointer = &newData;

  // copy all data from previous to new (to account for edges)
  for (size_t row = 0; row < newData.size(); ++row) {
    newData[row].resize(data[0].size());

    // each thread takes a row at a time
    for (size_t col = 0; col < newData[row].size(); ++col) {
      newData[row][col] = data[row][col];
    }
  }

  bool localEquilibrium = false;

  // while the matrix is not in equilibrium
  while (!localEquilibrium) {
    // run a stage and receive wrong cell amount
    wrongAmount = runStage(*dataPointer, *newDataPointer, jobInformation,
        jobInformation.equilibriumPointSentivity);

    // check if in equilibrium
    localEquilibrium = wrongAmount == 0;

    // if not in equilibrium
    if (!localEquilibrium) {
      /* swap matrixes (data will always have old data)
      * old data will be overwritten and lost within the matrixes
      */
      std::swap(*dataPointer, *newDataPointer);
    }

    // increase stages processed and reset coincidence amount
    stageCount++;
    wrongAmount = 0;
  }

  jobInformation.stateAmountRequired = stageCount;

  // write results on new file
  TermalTransferS::writeMatrixOnFile(newData, jobInformation);
}

// returns a matrix of the data found in the file data
int runStage(Matrix<double>& data, Matrix<double>& newData,
    JobInformation& jobInformation, const double equilibriumPointSentivity) {
  // get sizes
  size_t rowAmount = data.size(),
      colAmount = data[0].size();

  int notInEquilibrium = 0;

  // get formula parameters
  double time = jobInformation.stageTimeDuration;
  double alpha = jobInformation.termalDifusivity;
  double dimensions =
      jobInformation.cellDimensions * jobInformation.cellDimensions;

  // for each row
  #pragma omp for simd collapse (2)
  for (size_t row = 1; row < rowAmount - 1; ++row) {
    // for each column (each thread gets a row at a time)
    for (size_t col = 1; col < colAmount - 1; ++col) {
      // find the sum of neighbors
      double sumOfNeightbors =
          data[row][col - 1] + data[row - 1][col] +
          data[row][col + 1] + data[row + 1][col] -
          (4 * (data[row][col]));

      // find the new state of the block
      newData[row][col] = data[row][col] +
          (((time * alpha)/dimensions) * sumOfNeightbors);

      double difference = data[row][col] - newData[row][col];
      // ensure it is positive
      if (difference < 0) {
        difference = -difference;
      }

      // check if it is within epsilom
      if (difference > equilibriumPointSentivity) {
        notInEquilibrium++;
      }
    }
  }

  return notInEquilibrium;
}


// returns a matrix of the data found in the file data
int runStageST(Matrix<double>& data, Matrix<double>& newData,
JobInformation& jobInformation,
const double equilibriumPointSentivity) {
  // get sizes
  size_t rowAmount = data.size(),
      colAmount = data[0].size();

  int notInEquilibrium = 0;

  // get formula parameters
  double time = jobInformation.stageTimeDuration;
  double alpha = jobInformation.termalDifusivity;
  double dimensions =
      jobInformation.cellDimensions * jobInformation.cellDimensions;

  // for each row
  for (size_t row = 1; row < rowAmount - 1; ++row) {
    // for each column (each thread gets a row at a time)
    for (size_t col = 1; col < colAmount - 1; ++col) {
      // find the sum of neighbors
      double sumOfNeightbors =
          data[row][col - 1] + data[row - 1][col] +
          data[row][col + 1] + data[row + 1][col] -
          (4 * (data[row][col]));

      // find the new state of the block
      newData[row][col] = data[row][col] +
          (((time * alpha)/dimensions) * sumOfNeightbors);

      double difference = data[row][col] - newData[row][col];
      // ensure it is positive
      if (difference < 0) {
        difference = -difference;
      }

      // check if it is within epsilom
      if (difference > equilibriumPointSentivity) {
        notInEquilibrium++;
      }
    }
  }

  return notInEquilibrium;
}

// returns a matrix of the data found in the file data
void TermalTransferS::getMatrix
    (Matrix<double>& dataMatrix, JobInformation& jobInformation) {
  // get complete path with file name
  std::string fileName = jobInformation.extension + jobInformation.fileName;

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

// writes matrix data on binary file
void TermalTransferS::writeMatrixOnFile(Matrix<double>& dataMatrix,
    JobInformation& jobInformation) {
  // get dimensions
  size_t rowAmount = dataMatrix.size(),
      colAmount = dataMatrix[0].size();

  // get new file name
  std::string newFileName = jobInformation.extension +
      jobInformation.fileName.substr(0, jobInformation.fileName.size() - 4)
      + "-" + std::to_string(jobInformation.stateAmountRequired) +".bin";

  // create binary file
  std::ofstream newFile(newFileName, std::ios::binary);

  // ensure file is actually open
  if (!newFile.is_open()) {
    throw std::runtime_error
        ("Binary file could not be opened to write new plate");
  }

  // write matrix dimensions
  newFile.write(reinterpret_cast<const char*>(&rowAmount), sizeof(size_t));
  newFile.write(reinterpret_cast<const char*>(&colAmount), sizeof(size_t));

  // ensure size is correct
  dataMatrix.resize(rowAmount);

  // for each row
  for (size_t row = 0; row < rowAmount; ++row) {
    // write all columns
    newFile.write(reinterpret_cast<const char*>
        (&dataMatrix[row][0]), (sizeof(double) * colAmount));
  }

  std::cout << "File: " << newFileName << " written!" << std::endl;

  // close file
  newFile.close();
}

// erases all data related to JobInformation struct
void TermalTransferS::eraseJobData(std::vector<JobInformation>* jobData) {
  // erases job data
  delete(jobData);
}

// writes report of all simulations
void writeReport(std::vector<JobInformation>& jobsInformation,
    std::string fileName) {
  // creates report file name: e.g. filepath/job020.tsv
  std::string newFileName =
      fileName.substr(0, fileName.size() - 4) + ".tsv";

  // create binary file
  std::ofstream newFile(newFileName);

  size_t job = 0;
  // write all jobs into the report
  while (job < jobsInformation.size() &&
      newFile << jobsInformation[job++]) {
  }

  // close file
  newFile.close();
}
