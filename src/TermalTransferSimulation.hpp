#ifndef TERMALTRANSFERSIMULATION_HPP
#define TERMALTRANSFERSIMULATION_HPP

#include <iostream>
#include <vector>
#include <fstream>
#include <string>

struct JobInformation;


template <typename dataType>
using Matrix = std::vector<std::vector<dataType>>;

namespace TermalTransferS {
  std::vector<JobInformation>* getJobData(std::string& fileName);

  void processAllJobs(std::vector<JobInformation>* jobs);

  void processJob(JobInformation* jobInformation);

  Matrix<double>* getMatrix(JobInformation* jobInformation);

  void writeMatrixOnFile(Matrix<double>& dataMatrix, JobInformation* jobInformation);

  void eraseJobData(std::vector<JobInformation>*);
};

#endif