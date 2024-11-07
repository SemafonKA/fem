#pragma once
#include <string>
#include <sstream>
#include <fstream>

namespace fem::quads {

    //void fileReader(const std::string& filepath) {
    //    auto file = std::ifstream(filepath);

    //    if (!file.is_open()) {
    //        const auto errorStr = format("Error opening file \"{}\"\n", filepath);
    //        cerr << termcolor::red << errorStr << termcolor::reset;
    //        throw std::runtime_error(errorStr);
    //    }
    //    // Read entire file at once
    //    auto fileContent = std::string((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    //    file.close();
    //    auto iss = std::istringstream(fileContent);

    //    auto str = std::string();
    //    auto lineNumber = size_t(0);

    //    while (std::getline(iss, str)) {
    //        lineNumber++;
    //        auto charNumber = size_t(0);
    //    }
    //}

}
