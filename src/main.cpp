/* Copyright 2022 Joseph Stuart Valverde Kong.
Universidad de Costa Rica. CC BY 4.0 */
#include "TermalTransferSimulation.hpp"

int main(int argc, char* argv[]) {
  try {
    // run simulation
    TermalTransferS::runTermalTransferSimulation(argc, argv);
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << std::endl;
  }
}
