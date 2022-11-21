/* Copyright 2022 Joseph Stuart Valverde Kong.
Universidad de Costa Rica. CC BY 4.0 */
#include "TermalTransferSimulation.hpp"

int main(int argc, char* argv[]) {
  if (argc > 1) {
    // get file name
    std::string fileName(argv[1]);

    // get job data
    std::vector<JobInformation>* jobData =
    TermalTransferS::getJobData(fileName);

    // process everything
    TermalTransferS::processAllJobs(jobData, fileName);

    // erase everything
    TermalTransferS::eraseJobData(jobData);
  }
}
