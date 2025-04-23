#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept> // For runtime_error
#include <vector>
#include <map>       // For mapping column names to indices
#include <set>       // For efficient column checking
#include <algorithm> // For std::find, std::max_element
#include <limits>    // For numeric_limits

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
    int& num_samples)
{
    // --- 0. Initialization ---
    u.clear();
    y.clear();
    num_samples = 0;

    // --- 1. Open File ---
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filepath);
    }

    std::string line;
    std::vector<std::string> header;
    std::map<std::string, int> col_name_to_index;

    // --- 2. Read Header ---
    if (!std::getline(file, line)) {
        throw std::runtime_error("CSV file is empty or could not read header: " + filepath);
    }
    std::stringstream header_ss(line);
    std::string segment;
    int col_index = 0;
    while (std::getline(header_ss, segment, ',')) {
        // Trim whitespace (important for matching)
        segment.erase(0, segment.find_first_not_of(" \t\n\r\f\v"));
        segment.erase(segment.find_last_not_of(" \t\n\r\f\v") + 1);
        if (!segment.empty()) { // Avoid adding empty headers
           header.push_back(segment);
           col_name_to_index[segment] = col_index++;
        } else if (header_ss.peek() == ',') { // Handle consecutive commas creating empty columns
            // Treat as an unnamed column if necessary, or decide how to handle
             header.push_back(""); // Add empty string as placeholder if needed
             // col_name_to_index[""] = col_index++; // Or skip adding to map
             col_index++;
        }
    }
     if (header.empty()) {
         throw std::runtime_error("CSV header row is empty or invalid: " + filepath);
     }


    // --- 3. Validate and Map Requested Columns ---
    std::vector<int> u_col_indices; // Indices in the CSV corresponding to input_cols order
    std::vector<int> y_col_indices; // Indices in the CSV corresponding to output_cols order

    for (const auto& col_name : input_cols) {
        auto it = col_name_to_index.find(col_name);
        if (it == col_name_to_index.end()) {
            // Construct the error message listing available headers
            std::string available_headers = " Available headers are: ";
            for(size_t i = 0; i < header.size(); ++i) {
                available_headers += "'" + header[i] + "'" + (i == header.size() - 1 ? "." : ", ");
            }
            throw std::runtime_error("Requested input column '" + col_name + "' not found in CSV header." + available_headers);
        }
        u_col_indices.push_back(it->second);
    }

    for (const auto& col_name : output_cols) {
         auto it = col_name_to_index.find(col_name);
        if (it == col_name_to_index.end()) {
            // Construct the error message listing available headers
             std::string available_headers = " Available headers are: ";
            for(size_t i = 0; i < header.size(); ++i) {
                available_headers += "'" + header[i] + "'" + (i == header.size() - 1 ? "." : ", ");
            }
            throw std::runtime_error("Requested output column '" + col_name + "' not found in CSV header." + available_headers);
        }
        y_col_indices.push_back(it->second);
    }

    // --- 4. Read Data Rows into Temporary Row-Major Storage ---
    std::vector<std::vector<double>> u_rows; // Store rows temporarily
    std::vector<std::vector<double>> y_rows; // Store rows temporarily
    int line_num = 1; // Start counting after header

    while (std::getline(file, line)) {
        line_num++;
        if (line.empty() || line.find_first_not_of(" \t\n\r\f\v") == std::string::npos) {
            continue; // Skip empty lines
        }

        std::stringstream line_ss(line);
        std::vector<std::string> current_row_values;
        std::string value_str;

        // Read all values in the current row
        while (std::getline(line_ss, value_str, ',')) {
            current_row_values.push_back(value_str);
        }

        // Check if the row has enough columns for the requested indices
        int max_u_idx = u_col_indices.empty() ? -1 : *std::max_element(u_col_indices.begin(), u_col_indices.end());
        int max_y_idx = y_col_indices.empty() ? -1 : *std::max_element(y_col_indices.begin(), y_col_indices.end());
        int max_req_idx = std::max(max_u_idx, max_y_idx);

        if (max_req_idx >= 0 && (int)current_row_values.size() <= max_req_idx) {
             throw std::runtime_error("Error: Line " + std::to_string(line_num) + " in file " + filepath
                                      + " has only " + std::to_string(current_row_values.size())
                                      + " columns, but requested column index " + std::to_string(max_req_idx) + ".");
        }

        // Extract required values for U and Y for this row
        std::vector<double> current_u_row_doubles;
        std::vector<double> current_y_row_doubles;

        try {
            for (int csv_col_idx : u_col_indices) {
                std::string& val_str_ref = current_row_values[csv_col_idx];
                // Trim whitespace before conversion
                val_str_ref.erase(0, val_str_ref.find_first_not_of(" \t\n\r\f\v"));
                val_str_ref.erase(val_str_ref.find_last_not_of(" \t\n\r\f\v") + 1);
                current_u_row_doubles.push_back(std::stod(val_str_ref));
            }
            for (int csv_col_idx : y_col_indices) {
                 std::string& val_str_ref = current_row_values[csv_col_idx];
                 // Trim whitespace before conversion
                 val_str_ref.erase(0, val_str_ref.find_first_not_of(" \t\n\r\f\v"));
                 val_str_ref.erase(val_str_ref.find_last_not_of(" \t\n\r\f\v") + 1);
                current_y_row_doubles.push_back(std::stod(val_str_ref));
            }
        } catch (const std::invalid_argument& e) {
             throw std::runtime_error("Error: Non-numeric value encountered parsing line "
                                      + std::to_string(line_num) + " in file " + filepath + ". Details: " + e.what());
        } catch (const std::out_of_range& e) {
             throw std::runtime_error("Error: Numeric value out of range parsing line "
                                      + std::to_string(line_num) + " in file " + filepath + ". Details: " + e.what());
        }

        u_rows.push_back(current_u_row_doubles);
        y_rows.push_back(current_y_row_doubles);
        num_samples++;
    }

    file.close();

    // --- 5. Convert from Row Storage to Column-Major Vectors ---
    if (num_samples == 0) {
        std::cerr << "Warning: No valid data rows found in file: " << filepath << std::endl;
        return true; // No data, but technically successful read
    }

    size_t num_u_cols = input_cols.size();
    size_t num_y_cols = output_cols.size();

    // Resize final output vectors
    u.resize((size_t)num_samples * num_u_cols);
    y.resize((size_t)num_samples * num_y_cols);

    // Populate U (column-major)
    for (size_t j = 0; j < num_u_cols; ++j) { // Iterate through columns of U
        for (int i = 0; i < num_samples; ++i) { // Iterate through rows (samples)
            // Column-major index: row + col * num_rows
            u[i + j * num_samples] = u_rows[i][j];
        }
    }

    // Populate Y (column-major)
    for (size_t j = 0; j < num_y_cols; ++j) { // Iterate through columns of Y
        for (int i = 0; i < num_samples; ++i) { // Iterate through rows (samples)
            // Column-major index: row + col * num_rows
            y[i + j * num_samples] = y_rows[i][j];
        }
    }

    return true;
}

#endif // TEST_UTILS_H
