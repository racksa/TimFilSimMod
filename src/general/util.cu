// util.cpp

#include "util.hpp"
#include "../../globals.hpp"

bool hasEnding (std::string const &fullString, std::string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

std::unordered_map<std::string, std::unordered_map<std::string, std::string>> parseINI(const std::string& filename) {
    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> data;

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening the file: " << filename << std::endl;
        return data;
    }

    std::string current_section;
    std::string line;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key, value;

        if (line.empty() || line[0] == ';' || line[0] == '#') {
            // Skip empty lines and comments
            continue;
        } else if (line[0] == '[' && line[line.length() - 1] == ']') {
            // Section line (e.g., [SectionName])
            current_section = line.substr(1, line.length() - 2);
        } else if (std::getline(iss, key, '=') && std::getline(iss, value)) {
            // Key-value pair line
            data[current_section][key] = value;
        }
    }

    return data;
}

std::string data_from_ini(std::string section, std::string variable){
    std::string filename = "globals.ini";
    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> iniData = parseINI(filename);

    // Access and use the parameters
    if (iniData.find(section) != iniData.end()) {
        std::cout << section << " " << variable << ": " << iniData[section][variable] << std::endl;
        return iniData[section][variable];        
    } else {
        std::cout << "Section " << section << " not found in the INI file." << std::endl;
        return "";
    }
}



CartesianCoordinates spherical_to_cartesian(double r, double theta, double phi) {
    CartesianCoordinates cartesian;
    cartesian.x = r * sin(theta) * cos(phi);
    cartesian.y = r * sin(theta) * sin(phi);
    cartesian.z = r * cos(theta);
    return cartesian;
}