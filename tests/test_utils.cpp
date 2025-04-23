#include "test_utils.h"

bool load_test_data_from_csv(const std::string& filepath, 
                           std::vector<double>& u, 
                           std::vector<double>& y) {
    // Open the CSV file
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filepath << std::endl;
        return false;
    }
    
    // Clear output vectors
    u.clear();
    y.clear();
    
    // Read the header line
    std::string line;
    if (!std::getline(file, line)) {
        std::cerr << "Empty file: " << filepath << std::endl;
        return false;
    }
    
    // Read data lines
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string u_value, y_value;
        
        // Parse the U value
        if (!std::getline(ss, u_value, ',')) {
            continue;
        }
        
        // Parse the Y value
        if (!std::getline(ss, y_value, ',')) {
            continue;
        }
        
        try {
            // Convert and store values
            u.push_back(std::stod(u_value));
            y.push_back(std::stod(y_value));
        } catch (const std::exception& e) {
            std::cerr << "Error parsing line: " << line << " - " << e.what() << std::endl;
            // Continue processing other lines
        }
    }
    
    // Check if we got data
    if (u.empty() || y.empty()) {
        std::cerr << "No valid data found in file: " << filepath << std::endl;
        return false;
    }
    
    return true;
}
