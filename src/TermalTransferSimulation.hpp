/* Copyright 2022 Joseph Stuart Valverde Kong.
Universidad de Costa Rica. CC BY 4.0 */
#ifndef TERMALTRANSFERSIMULATION_HPP
#define TERMALTRANSFERSIMULATION_HPP

#include <iostream>
#include <vector>
#include <fstream>
#include <string>

/**
 * @brief forward declaration of JobInformation
 * 
 */
struct JobInformation;

/**
 * @brief defines a vector with vectors of a type as a matrix of a type
 * 
 * @tparam dataType 
 */
template <typename dataType>
using Matrix = std::vector<std::vector<dataType>>;

namespace TermalTransferS {
  /**
   * @brief returns a vector of jobs found within file from file name
   * 
   * @param fileName name of file where jobs are extracted
   * @return std::vector<JobInformation>* 
   */
  std::vector<JobInformation>* getJobData(std::string& fileName);

  /**
   * @brief processess all jobs in the vector
   * 
   * @param jobs vector of jobs to be processed
   * @param fileName name of file where jobs were extracted
   */
  void processAllJobs(std::vector<JobInformation>* jobs, std::string& fileName);

  /**
   * @brief processess the given job
   * 
   * @param jobInformation job to be processed
   */
  void processJob(JobInformation* jobInformation);

  /**
   * @brief returns a matrix of the data found in the file data
   * @details opens file from jobInformation filename to extract data
   * @param jobInformation has file name of file with binary data
   * @return Matrix<double>* 
   */
  Matrix<double>* getMatrix(JobInformation* jobInformation);

  /**
   * @brief writes matrix data on binary file
   * 
   * @param dataMatrix matrix with data to be written
   * @param jobInformation information for file name and 
   * other related information
   */
  void writeMatrixOnFile(Matrix<double>& dataMatrix,
  JobInformation* jobInformation);

  /**
   * @brief erases all data related to JobInformation struct
   * 
   */
  void eraseJobData(std::vector<JobInformation>*);
};  // namespace TermalTransferS

#endif
