/* Copyright 2022 Joseph Stuart Valverde Kong.
Universidad de Costa Rica. CC BY 4.0 */
#ifndef TERMALDATA_HPP
#define TERMALDATA_HPP

#include <string>
#include <fstream>

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

  /**
   * @brief overload of operator >>
   * @details reads from a file all information for a job
   * @param file from where the data will be read
   * @param jobInformation where the data will be stored
   * @return std::ifstream& 
   */
  friend std::ifstream& operator>> (std::ifstream& file,
      JobInformation& jobInformation) {
    // read from file and into jobInformation all data
    file >> jobInformation.fileName;
    file >> jobInformation.stageTimeDuration;
    file >> jobInformation.termalDifusivity;
    file >> jobInformation.cellDimensions;
    file >> jobInformation.equilibriumPointSentivity;

    return file;
  }

  /**
   * @brief overload of write operator
   * @details writes all information from a job into the report file
   * @param file where the report will be written
   * @param jobInformation from where the information is extracted
   * @return std::ofstream& 
   */
  friend std::ofstream& operator<< (std::ofstream& file,
      JobInformation& jobInformation) {
    // calculate total time
    double totalTime = jobInformation.stageTimeDuration
        * jobInformation.stateAmountRequired;

    // write all information into file
    file << jobInformation.fileName
        << "\t" << jobInformation.stageTimeDuration
        << "\t" << jobInformation.termalDifusivity
        << "\t" << jobInformation.cellDimensions
        << "\t" << jobInformation.equilibriumPointSentivity
        << "\t" << jobInformation.stateAmountRequired
        << "\t" << format_time(totalTime) << "\n";

    return file;
  }

 private:
  /**
   * @brief gets the time formated as string
   * 
   * @param seconds amount of seconds taken
   * @return std::string
   */
  static std::string format_time(const time_t seconds) {
    char text[48];  // YYYY/MM/DD hh:mm:ss
    const std::tm& gmt = * std::gmtime(&seconds);
        snprintf(text, sizeof(text),
        "%04d/%02d/%02d\t%02d:%02d:%02d", gmt.tm_year
        - 70, gmt.tm_mon, gmt.tm_mday - 1, gmt.tm_hour,
        gmt.tm_min, gmt.tm_sec);
    return text;
  }
};

#endif
