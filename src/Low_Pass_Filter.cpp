// Copyright (C) 2018, Keith Pedersen (Keith.David.Pedersen@gmail.com)
// All rights reserved

// This file calls h_Gaussian for the supplied lambda angles and writes hl to a file.

// CM: Edits by Mithila:
                        // - changing the ini to conf
			// - adding output file handling and variable lMax
			// - changing how Qt hadnled the string of lambda values (the original code's implementation of this failed with current libraries)
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include <QtCore/QSettings>
#include <QtCore/QStringList>

#include <sys/stat.h>

#include "kdpTools.hpp"
#include "ShapeFunction.hpp"

// ------------------------------------------------------------
// Helper: check if directory exists
// ------------------------------------------------------------
static bool DirExists(const std::string& path)
{
    struct stat info;
    return (stat(path.c_str(), &info) == 0 &&
            (info.st_mode & S_IFDIR));
}

// ------------------------------------------------------------
// Main
// ------------------------------------------------------------
int main()
{
    QSettings parsedINI("/home/mithila/EmPowerSpectrum/CONF/Low_Pass_Filter.conf", QSettings::IniFormat);
    parsedINI.sync();
    std::cerr << "Reading config file: "
	      << parsedINI.fileName().toStdString() << "\n";
    

    // --------------------------------------------------------
    // Output directory
    // --------------------------------------------------------
    std::string outputDir =
        parsedINI.value("main/outputDir", "./").toString().toStdString();

    if (!outputDir.empty() &&
        outputDir.back() != '/' &&
        outputDir.back() != '\\')
    {
        outputDir += "/";
    }

    if (!DirExists(outputDir)) {
        std::cerr << "Warning: output directory does not exist: "
                  << outputDir
                  << "\nUsing current directory instead.\n";
        outputDir = "./";
    }

    // --------------------------------------------------------
    // Read lambda list (comma-separated)
    // --------------------------------------------------------
    QStringList lambdaValues =
      parsedINI.value("main/lambda").toStringList();
    
    for (auto& s : lambdaValues)
      s = s.trimmed();
    
    if (lambdaValues.empty()) {
      std::cerr << "Error: no lambda values provided\n";
      return 1;
    }
    
    // --------------------------------------------------------
    // lMax
    // --------------------------------------------------------
    size_t const lMax =
        size_t(parsedINI.value("main/lMax", 512).toUInt());

    // --------------------------------------------------------
    // Build shape functions
    // --------------------------------------------------------
    std::vector<h_Gaussian> shapeVec;
    shapeVec.reserve(lambdaValues.size());

    for (auto const& lambdaQString : lambdaValues) {
        double const lambda =
            kdp::ReadAngle<double>(lambdaQString.toStdString());
        shapeVec.emplace_back(lambda);
    }

    // --------------------------------------------------------
    // Write output file
    // --------------------------------------------------------
    std::string outputFile = outputDir + "hl.dat";
    std::ofstream file(outputFile, std::ios::trunc);

    if (!file.is_open()) {
        std::cerr << "Error: could not open output file "
                  << outputFile << std::endl;
        return 1;
    }

    char buff[1024];

    // Header
    file << "#";
    for (auto const& lambdaQString : lambdaValues)
        file << "  " << lambdaQString.toStdString();
    file << "\n";

    // l = 0 line
    sprintf(buff, "%4lu", 0lu);
    file << buff;

    for (size_t i = 0; i < shapeVec.size(); ++i) {
        sprintf(buff, "  %.16e", 1.0);
        file << buff;
    }
    file << "\n";

    // l >= 1
    for (size_t l = 1; l <= lMax; ++l) {
        sprintf(buff, "%4lu", l);
        file << buff;

        for (auto& shape : shapeVec) {
            sprintf(buff, "  %.16e", shape.hl(l));
            file << buff;
        }
        file << "\n";
    }

    std::cout << "Low-pass filter written to:\n  "
              << outputFile << std::endl;

    return 0;
}
