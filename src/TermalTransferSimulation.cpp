/* Copyright 2022 Joseph Stuart Valverde Kong.
Universidad de Costa Rica. CC BY 4.0 */
#include <cstdlib>
#include <utility>
#include <ctime>
#include <iomanip>
#include <omp.h>
#include <thread>
#include <mpi.h>

#include "TermalTransferSimulation.hpp"
#include "TermalData.hpp"

/**
 * @brief administrative process handling task issue and process coordination
 * 
 * @param fileName 
 * @param size 
 */
void mpiRank0(std::string fileName, const int32_t size);

/**
 * @brief sends data to all processes that request it
 * 
 * @param jobData vector of all jobs
 * @param stringSizes sizes of the strings to be sent
 * @param dataToSend 
 * @param jobLocationArr 
 */
void sendData(std::vector<JobInformation>& jobData, int32_t* stringSizes,
    double* dataToSend, int32_t* jobLocationArr);

/**
 * @brief sends stop conditions to stop all processes
 * @details stop condition is fileName size 0
 * @param stringSizes sizes of strings
 * @param size amount of processes
 */
void stopProcessess(int32_t* stringSizes, const int32_t size);

/**
 * @brief worker processes assigned to processing
 * 
 * @param rank process number or ID
 * @param threadAmount amount of desired threads
 */
void mpiRankAny(const int32_t rank, const int32_t threadAmount);

/**
 * @brief Receives data for worker process from administrating process
 * 
 * @param jobs vector where jobs will be stored
 * @param processedJobs amount of jobs already processed
 * @param stringSizes size of strings to be received
 * @param jobData buffer where job data is stored before being sent to the job
 * @param rank process number or ID
 * @return int32_t 
 */
int32_t receiveData(std::vector<JobInformation>& jobs, int32_t& processedJobs,
    int32_t* stringSizes, double* jobData, const int32_t rank);


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
JobInformation* jobInformation,
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
      return;
    }

    int rank = 0, size = 0;

    // get the rank of the current process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // get the amount of processes
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // if there is just one process
    if (size == 1) {
      // get job data
      std::vector<JobInformation>* jobData =
      TermalTransferS::getJobData(fileName);

      // process everything
      TermalTransferS::processAllJobs(jobData, fileName, threadAmount);

      // erase everything
      TermalTransferS::eraseJobData(jobData);

    // if there are multiple process
    } else {
      // if rank 0
      if (rank == 0) {
        // (administrating process)
        mpiRank0(fileName, size);
      // for all other process
      } else {
        // (working process)
        mpiRankAny(rank, threadAmount);
      }
    }

    // finalize all process
    MPI_Finalize();
  }
}

void mpiRank0(std::string fileName, const int32_t size) {
  // get jobs
  std::vector<JobInformation>* jobData =
      TermalTransferS::getJobData(fileName);
  // allocate space for string size vector
  int32_t* stringSizes = (int32_t*) calloc(2, sizeof(int32_t));
  // allocate space for the data parameters of each job
  double* dataToSend = (double*) calloc(4, sizeof(double));
  // allocate space for array holding which job is being processed by which process
  int32_t* jobLocationArr = (int32_t*) calloc(jobData->size(), sizeof(int32_t));

  // send data to all process
  sendData(*jobData, stringSizes, dataToSend, jobLocationArr);

  // stop processess with stop conditions
  stopProcessess(stringSizes, size);

  // tell other processes to send times
  for (int process = 1; process < size; ++process) {
    MPI_Send(&process, /* amount */1, MPI_INT,
        process, 0, MPI_COMM_WORLD);
  }

  // for all jobs (receive times)
  for (size_t job = 0; job < jobData->size(); ++job) {
    // receive time
    MPI_Recv(&(*jobData)[job].stateAmountRequired, 1,
    MPI_INT, jobLocationArr[job], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  // write report
  writeReport(*jobData, fileName);

  // erase data
  TermalTransferS::eraseJobData(jobData);
  free(jobLocationArr);
}

void sendData(std::vector<JobInformation>& jobData, int32_t* stringSizes,
    double* dataToSend, int32_t* jobLocationArr) {
  int32_t rankToSend = -1;
  // get size of extension, it is shared among all
  stringSizes[1] = jobData[0].extension.size();

  // for all jobs
  for (size_t job = 0; job < jobData.size(); ++job) {
    // receive signal from processes ready to process
    MPI_Recv(&rankToSend, 1,
    MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    jobLocationArr[job] = rankToSend;

    // send sizes of name and extension
    stringSizes[0] = jobData[job].fileName.size();
    MPI_Send(stringSizes, /* amount */2, MPI_INT,
        rankToSend, 0, MPI_COMM_WORLD);

    // send name
    MPI_Send(&jobData[job].fileName.c_str()[0], /* amount */stringSizes[0], MPI_UNSIGNED_CHAR,
        rankToSend, 0, MPI_COMM_WORLD);

    // send extension
    MPI_Send(&(jobData[job].extension.c_str()[0]), /* amount */stringSizes[1],
     MPI_UNSIGNED_CHAR, rankToSend, 0, MPI_COMM_WORLD);

    // send 4 datas
    dataToSend[0] = jobData[job].stageTimeDuration;
    dataToSend[1] = jobData[job].termalDifusivity;
    dataToSend[2] = jobData[job].cellDimensions;
    dataToSend[3] = jobData[job].equilibriumPointSentivity;

    MPI_Send(dataToSend, /* amount */4, MPI_DOUBLE,
        rankToSend, 0, MPI_COMM_WORLD);
  }
}

void stopProcessess(int32_t* stringSizes, const int32_t size) {
  int32_t rankToSend = -1;

  // for all processess (stop conditions)
  for (int process = 0; process < size - 1; ++process) {
    // receive ready to process
    MPI_Recv(&rankToSend, 1,
    MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // signal there is no more to do (string sizes = 0)
    stringSizes[0] = 0;
    stringSizes[1] = 0;

    MPI_Send(&stringSizes[0], /* amount */2, MPI_INT,
        rankToSend, 0, MPI_COMM_WORLD);
  }
}

void mpiRankAny(const int32_t rank, const int32_t threadAmount) {
  std::vector<JobInformation> jobs;

  int32_t processedJobs = 0;
  // buffer to sizes of file name and extension to be received
  int32_t* stringSizes = (int32_t*) calloc(2, sizeof(int32_t));
  // buffer for received job data
  double* jobData = (double*) calloc(4, sizeof(double));

  // while there are jobs
  while(true) {
    // receive all data for the job
    int dataState =
    receiveData(jobs, processedJobs, stringSizes, jobData, rank);

    // end loop if no more data has been received
    if (dataState == 1) {
      break;
    }

    // process job
    TermalTransferS::processJob(&jobs[processedJobs], threadAmount);

    processedJobs++;
  }

  // wait until told to continue
  int buffer = 0;
  MPI_Recv(&buffer, 1,
    MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  // std::cout << rank << " done processing: " << processedJobs << std::endl;
  // for all processed jobs
  for (int job = 0; job < processedJobs; ++job) {
    // send the time
    MPI_Send(&jobs[job].stateAmountRequired, /* amount */1, MPI_INT,
        0, 0, MPI_COMM_WORLD);
  }

  free(stringSizes);
  free(jobData);
}

int32_t receiveData(std::vector<JobInformation>& jobs, int32_t& processedJobs,
    int32_t* stringSizes, double* jobData, const int32_t rank) {

  // send ready to process
  MPI_Send(&rank, 1, MPI_INT, /*destination*/ 0, 0, MPI_COMM_WORLD);

  // receive string sizes
  MPI_Recv(stringSizes, 2,
    MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  // if file name size == 0
  if (stringSizes[0] == 0) {
    // break cycle, no more to do
    return 1;
  }
  // add new job to the vector
  jobs.emplace_back(JobInformation());

  // receive fileName
  std::string fileName(stringSizes[0], '\0');
  MPI_Recv(&fileName[0], stringSizes[0],
    MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  
  jobs[processedJobs].fileName = fileName;

  // receive extension
  std::string extension(stringSizes[1], '\0');
  MPI_Recv(&extension[0], stringSizes[1],
    MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  jobs[processedJobs].extension = extension;

  // receive data
  MPI_Recv(jobData, /* amount */4, MPI_DOUBLE,
      0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  jobs[processedJobs].stageTimeDuration = jobData[0];
  jobs[processedJobs].termalDifusivity = jobData[1];
  jobs[processedJobs].cellDimensions = jobData[2];
  jobs[processedJobs].equilibriumPointSentivity = jobData[3];

  return 0;
}

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

    // place the same extension on all JobInformations
    (*dataVector)[jobsFoundAmount].extension = extension;

    dataPosition++;

    // advance to next job if all current data already placed
    if (dataPosition == 5) {
      jobsFoundAmount++;
      dataPosition = 0;
    }
  }

  // shrink to size
  dataVector->resize(jobsFoundAmount);

  return dataVector;
}

// processess all jobs in the vector
void TermalTransferS::processAllJobs(std::vector<JobInformation>* jobs,
    std::string& fileName, int32_t threadAmount) {
  int jobsAmount = jobs->size();

  // for all jobs
  for (int currentJob = 0; currentJob < jobsAmount; ++currentJob) {
    // process them
    processJob(&(*jobs)[currentJob], threadAmount);
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
void TermalTransferS::processJob
    (JobInformation* jobInformation, int32_t threadAmount) {
  // get the matrix out of the file
  Matrix<double>* data = getMatrix(jobInformation);

  // set a new matrix for new data
  Matrix<double>* newData = new Matrix<double>;

  newData->resize(data->size());

  size_t stageCount = 0;

  bool inEquilibrium = false;

  size_t wrongAmount = 0;

  #pragma omp parallel num_threads(threadAmount) \
    default(none) shared(data, newData, jobInformation, threadAmount,\
    stageCount, inEquilibrium, wrongAmount)
  {
    // copy all data from previous to new (to account for edges)
    #pragma omp for
    for (size_t row = 0; row < newData->size(); ++row) {
      (*newData)[row].resize((*data)[0].size());

      // each thread takes a row at a time
      for (size_t col = 0; col < (*newData)[row].size(); ++col) {
        (*newData)[row][col] = (*data)[row][col];
      }
    }

    bool localEquilibrium = false;

    // while the matrix is not in equilibrium
    while (!localEquilibrium) {
      // run a stage and receive wrong cell amount
      int localWrong = runStage(*data, *newData, jobInformation,
      jobInformation->equilibriumPointSentivity);

      #pragma omp critical
      {
        wrongAmount += localWrong;
      }

      #pragma omp barrier

      // check if in equilibrium
      #pragma omp single
      inEquilibrium = wrongAmount == 0;

      #pragma omp critical
      {
        localEquilibrium = inEquilibrium;
      }
      
      // if not in equilibrium
      #pragma omp single
      if (!localEquilibrium) {
        /* swap matrixes (data will always have old data)
        * old data will be overwritten and lost within the matrixes
        */
        std::swap(data, newData);
      }

      // increase stages processed and reset coincidence amount
      #pragma omp single
      {
        stageCount++;
        wrongAmount = 0;
      }
    }
  }

  jobInformation->stateAmountRequired = stageCount;

  // write results on new file
  writeMatrixOnFile(*newData, jobInformation);
  delete(data);
}

// returns a matrix of the data found in the file data
int runStage(Matrix<double>& data, Matrix<double>& newData,
JobInformation* jobInformation,
const double equilibriumPointSentivity) {
  // get sizes
  size_t rowAmount = data.size(),
      colAmount = data[0].size();

  int notInEquilibrium = 0;

  // get formula parameters
  double time = jobInformation->stageTimeDuration,
      alpha = jobInformation->termalDifusivity,
      dimensions =
      jobInformation->cellDimensions * jobInformation->cellDimensions;

  // for each row
  #pragma omp for
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
Matrix<double>* TermalTransferS::getMatrix(JobInformation* jobInformation) {
  // get complete path with file name
  std::string fileName = jobInformation->extension + jobInformation->fileName;

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
