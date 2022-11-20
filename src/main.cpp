#include "TermalTransferSimulation.hpp"

int main(int argc, char* argv[]) {
  if (argc > 1) {
    std::string fileName(argv[1]);

    std::vector<JobInformation*>* jobData =
    TermalTransferS::getJobData(fileName);

    TermalTransferS::processAllJobs(jobData);
  }
}