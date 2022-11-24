/* Copyright 2022 Joseph Stuart Valverde Kong.
Universidad de Costa Rica. CC BY 4.0 */
#ifndef TERMALDATA_HPP
#define TERMALDATA_HPP

#include <string>

struct JobInformation {
  /**
   * @brief name of the binary file
   */
  std::string fileName;

  /**
   * @brief extension to lead to binary file
   */
  std::string extension;

  /**
   * @brief how much time each stage is to last
   */
  double stageTimeDuration = 0;

  /**
   * @brief terman difusivity of material where heat is spread
   */
  double termalDifusivity = 0;

  /**
   * @brief dimensions of each cell holding heat
   */
  double cellDimensions = 0;

  /**
   * @brief epsilon limit where anything under
   * is considered as under equilibrium
   */
  double equilibriumPointSentivity = 0;

  /**
   * @brief total amount of stages that it took for plate
   * to reach equilibrium
   */
  int stateAmountRequired = 1;

  JobInformation(){}
};

#endif
