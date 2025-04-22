#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <filesystem>

class CSVReader {
private:
    std::string filePath;
    std::vector<std::string> headers;
    std::vector<std::vector<std::string>> data;
    bool isLoaded = false;

    

    // Helper function to convert string to type T
    template<typename T>
    T convert(const std::string& value) const {
        std::istringstream iss(value);
        T result;
        if (!(iss >> result)) {
            throw std::runtime_error("Failed to convert value: " + value);
        }
        return result;
    }

    // Specialization for std::string (no conversion needed)
    template<>
    std::string convert<std::string>(const std::string& value) const {
        return value;
    }
    
    void loadData() 
    {
        std::ifstream file(filePath);
        if (!file.is_open()) {
            throw std::runtime_error("Error opening file at: " + filePath);
        }

        std::string line;
        // Read headers
        if (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string value;
            while (std::getline(ss, value, ',')) {
                headers.push_back(value);
            }
        }

        // Read data rows
        while (std::getline(file, line)) {
            std::vector<std::string> row;
            std::stringstream ss(line);
            std::string value;
            while (std::getline(ss, value, ',')) {
                row.push_back(value);
            }
            data.push_back(row);
        }

        isLoaded = true;
    }

public:

    // Constructor with hardcoded Mac path
    CSVReader() : filePath("/Users/hin601/Documents/TestBuild/CRHM_Compare/PythonScripts/unittestscripts/Cleaned_data.csv") {
        if (!std::filesystem::exists(filePath)) {
            throw std::runtime_error("File not found: " + filePath);
        }

        loadData();
    }

    
    
    // Function to get value from CSV
    template<typename T>
    T getValue(const std::string& columnName, int rowNumber) 
    {
        if (!isLoaded)
        {
            throw std::runtime_error("Data not loaded");
        }
        
        // Find the target column
        auto it = std::find(headers.begin(), headers.end(), columnName);
        if (it == headers.end()) 
        {
            throw std::runtime_error("Column '" + columnName + "' not found.");
        }
        int targetColumn = std::distance(headers.begin(), it);
        
        // Check row bounds
        if (rowNumber < 0 || rowNumber >= static_cast<int>(data.size())) 
        {
            throw std::runtime_error("Row " + std::to_string(rowNumber) + " is out of bounds");
        }

        // Check column bounds
        if (targetColumn >= static_cast<int>(data[rowNumber].size())) 
        {
            throw std::runtime_error("Column index out of bounds for row " + 
                                   std::to_string(rowNumber));
        }

        return convert<T>(data[rowNumber][targetColumn]);
    } 
        
        //std::ifstream file(filePath);
        //std::string line, value;
        //std::vector<std::string> columnNames;
        //int targetColumn = -1;
        //int currentRow = 0;

        //// Check if file opened successfully
        //if (!file.is_open()) {
        //    throw std::runtime_error("Error opening file at: " + filePath);
        //}

        //// Read header row to get column names
        //if (std::getline(file, line)) {
        //    std::stringstream ss(line);
        //    while (std::getline(ss, value, ',')) {
        //        columnNames.push_back(value);
        //    }

        //    // Find the target column
        //    for (size_t i = 0; i < columnNames.size(); ++i) {
        //        if (columnNames[i] == columnName) {
        //            targetColumn = static_cast<int>(i);
        //            break;
        //        }
        //    }

        //    if (targetColumn == -1) {
        //        throw std::runtime_error("Column '" + columnName + "' not found.");
        //    }
        //}

        //// Read data rows until we reach the target row
        //while (std::getline(file, line)) {
        //    if (currentRow == rowNumber) {
        //        std::vector<std::string> rowData;
        //        std::stringstream ss(line);
        //        while (std::getline(ss, value, ',')) {
        //            rowData.push_back(value);
        //        }

        //        if (targetColumn < rowData.size()) {
        //            return convert<T>(rowData[targetColumn]);
        //        } else {
        //            throw std::runtime_error("Column index out of bounds for row " + std::to_string(rowNumber));
        //        }
        //    }
        //    currentRow++;
        //}
        //
        //throw std::runtime_error("Row " + std::to_string(rowNumber) + " not found in the CSV file.");
    //}
};
