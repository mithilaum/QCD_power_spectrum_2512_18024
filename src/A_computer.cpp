#include <QtCore/QSettings>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

#include "PowerSpectrum.hpp"
#include "ShapeFunction.hpp"
#include "Admin.hpp"
#include "kdpVectors.hpp"

// ------------------------------------------------------------
// Helpers
// ------------------------------------------------------------

static std::string CoeffToString(double x)
{
    std::ostringstream ss;
    ss << x;   // NO fixed precision, as requested
    return ss.str();
}

// ------------------------------------------------------------
// Read final-state particles from PYTHIA .dat file
// ------------------------------------------------------------

static std::vector<PowerSpectrum::PhatF>
ReadFinalStates(const std::string& filePath)
{
    std::ifstream file(filePath);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open input data file: " + filePath);
    }

    std::vector<PowerSpectrum::PhatF> particles;
    std::string line;

    // skip header
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::istringstream iss(line);

        int idx, pdg, status;
        double E, px, py, pz, m;
        int charge, visible;

        if (!(iss >> idx >> pdg >> status >> E >> px >> py >> pz >> m >> charge >> visible))
            continue;

        if (!visible)
            continue;

        particles.emplace_back(px, py, pz, E);
    }

    return PowerSpectrum::PhatF::To_PhatF_Vec(particles);
}

// ------------------------------------------------------------
// Main
// ------------------------------------------------------------

int main(int argc, char** argv)
{
    const std::string confPath =
        (argc > 1) ? argv[1] : "./CONF/Angular_Correlation.conf";

    QSettings parsedINI(confPath.c_str(), QSettings::IniFormat);

    // -----------------------------
    // Read configuration
    // -----------------------------

    std::string inputDataFile =
        parsedINI.value("main/inputDataFile", "").toString().toStdString();

    if (inputDataFile.empty()) {
        std::cerr << "Error: inputDataFile not specified\n";
        return 1;
    }

    std::string outputDir =
        parsedINI.value("main/outputDir", "./").toString().toStdString();

    size_t zSamples =
        parsedINI.value("main/zSamples", 2048).toUInt();

    size_t lMax =
        parsedINI.value("power/lMax", 1024).toUInt();

    double underCoeff =
        parsedINI.value("smear/Coeff_under_smeared", 0.1).toDouble();

    double overCoeff =
        parsedINI.value("smear/Coeff_over_smeared", 10.0).toDouble();

    double trackRes =
        kdp::ReadAngle<double>(
            parsedINI.value("smear/tracks", "1 deg").toString().toStdString()
        );

    // -----------------------------
    // Prepare output filename
    // -----------------------------

    if (!outputDir.empty() &&
        outputDir.back() != '/' &&
        outputDir.back() != '\\')
    {
        outputDir += "/";
    }

    if (!std::ifstream(outputDir).good()) {
        std::cerr << "Warning: output directory does not exist, using ./\n";
        outputDir = "./";
    }

    std::string inputName = GetProcess(inputDataFile);

    std::string outputFile =
        outputDir +
        "A_" + inputName + "_" +
        CoeffToString(underCoeff) + "_" +
        CoeffToString(overCoeff) +
        ".dat";

    // -----------------------------
    // Read particles
    // -----------------------------

    auto particles = ReadFinalStates(inputDataFile);

    // -----------------------------
    // Shape functions
    // -----------------------------

    h_Gaussian naturalShape(trackRes);
    h_Gaussian underShape(underCoeff * naturalShape.Lambda());
    h_Gaussian overShape(overCoeff * naturalShape.Lambda());

    // -----------------------------
    // Power spectrum
    // -----------------------------

    size_t numThreads = parsedINI.value("main/numThreads", 4).toUInt();
    PowerSpectrum H(numThreads);

    auto Hl_raw     = H.Hl_Obs(lMax, particles);
    auto Hl_nat     = H.Hl_Obs(lMax,
                          PowerSpectrum::ShapedParticleContainer(particles, naturalShape));
    auto Hl_under   = H.Hl_Obs(lMax,
                          PowerSpectrum::ShapedParticleContainer(particles, underShape));
    auto Hl_over    = H.Hl_Obs(lMax,
                          PowerSpectrum::ShapedParticleContainer(particles, overShape));

    // -----------------------------
    // Angular correlations
    // -----------------------------

    auto AC = PowerSpectrum::AngularCorrelation(
        { Hl_raw, Hl_nat, Hl_under, Hl_over },
        zSamples
    );

    auto const& zVals = AC.first;
    auto const& Avals = AC.second;

    // -----------------------------
    // Write output
    // -----------------------------

    std::ofstream out(outputFile);
    if (!out.is_open()) {
      std::cerr << "Error: could not open output file "
		<< outputFile << std::endl;
      return 1;
    }
    
    // Header
    out << "#                     z                      raw                 natural                    under                     over\n";

    // Match legacy numeric formatting
    out.setf(std::ios::scientific);
    out.precision(16);

    for (size_t i = 0; i < zVals.size(); ++i) {
      out << zVals[i];
      for (size_t j = 0; j < Avals.size(); ++j)
        out << "  " << Avals[j][i];
      out << "\n";
    }
    
    std::cout << "Angular correlation written to:\n  "
              << outputFile << std::endl;

    return 0;
}
