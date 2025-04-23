#ifndef SLICOT_TEST_UTILS_H
#define SLICOT_TEST_UTILS_H

#include <vector>
#include <string>
#include <stdexcept> // Include for std::runtime_error documentation

/**
 * @brief Loads specific columns from a CSV file into column-major vectors.
 *
 * Reads a CSV file, validates requested input and output column names against
 * the header, and loads the corresponding data into separate column-major vectors.
 * Assumes the first row of the CSV is the header and data is comma-separated.
 *
 * @param filepath Path to the CSV file.
 * @param input_cols Vector of strings containing the header names of columns to load into u.
 * @param output_cols Vector of strings containing the header names of columns to load into y.
 * @param u Output vector (will be resized and populated in column-major format) for input data.
 * Dimensions will be num_samples x input_cols.size().
 * @param y Output vector (will be resized and populated in column-major format) for output data.
 * Dimensions will be num_samples x output_cols.size().
 * @param num_samples Output parameter passed by reference, will store the number of data rows read.
 * @return True if loading was successful (file opened, headers matched, data parsed).
 * Note: Errors during processing typically throw exceptions rather than returning false.
 * @throws std::runtime_error If the file cannot be opened, a requested column header is not found,
 * a row has an unexpected number of columns, or if non-numeric data is encountered
 * in a requested column.
 */
bool load_test_data_from_csv(
    const std::string& filepath,
    const std::vector<std::string>& input_cols,
    const std::vector<std::string>& output_cols,
    std::vector<double>& u,
    std::vector<double>& y,
    int& num_samples);

#endif // SLICOT_TEST_UTILS_H
