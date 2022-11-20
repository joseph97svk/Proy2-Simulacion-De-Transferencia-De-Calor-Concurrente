#ifndef TERMALDATA_HPP
#define TERMALDATA_HPP

#include <string>

struct JobInformation {
  std::string fileName;
  std::string extension;
  double stageTimeDuration = 0;
  double termalDifusivity = 0;
  double cellDimensions = 0;
  double equilibriumPointSentivity = 0;

  JobInformation(){}
};

#endif