#ifndef SLICOT_TEST_UTILS_H
#define SLICOT_TEST_UTILS_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

/**
 * @brief Load test data from a CSV file with columns U and Y
 * 
 * @param filepath Path to the CSV file
 * @param u Vector to store U data (inputs)
 * @param y Vector to store Y data (outputs)
 * @return true if successful, false otherwise
 */
bool load_test_data_from_csv(const std::string& filepath, 
                           std::vector<double>& u, 
                           std::vector<double>& y);

#endif // SLICOT_TEST_UTILS_H
