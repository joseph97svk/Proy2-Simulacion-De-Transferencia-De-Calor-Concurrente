/* Copyright 2022 Joseph Stuart Valverde Kong.
Universidad de Costa Rica. CC BY 4.0 */
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
  int plateRepeatIndex = 1;

  JobInformation(){}
};

#endif
