#include "Post_PYTHIA_PowerJets.hpp"
#include "Admin.hpp"
#include "PowerSpectrum.hpp"
#include "ArrogantDetector.hpp"
#include "kdpVectors.hpp"
#include "kdpStdVectorMath.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <QtCore/QSettings>

// Calculates the angular separation between two 3-vectors using their dot product
// Returns the angle in radians between 0 and π
double angularDistance(const kdp::Vec3& v1, const kdp::Vec3& v2) {
    // Compute dot product of the vectors
    double dot = v1.x1 * v2.x1 + v1.x2 * v2.x2 + v1.x3 * v2.x3;
    
    // Calculate magnitudes of both vectors
    double mag1 = std::sqrt(v1.x1 * v1.x1 + v1.x2 * v1.x2 + v1.x3 * v1.x3);
    double mag2 = std::sqrt(v2.x1 * v2.x1 + v2.x2 * v2.x2 + v2.x3 * v2.x3);
    
    // Calculate cos(theta) = dot/(|v1||v2|) with numerical safety checks
    double cosTheta = dot / (mag1 * mag2);
    if (cosTheta > 1.0) cosTheta = 1.0;
    if (cosTheta < -1.0) cosTheta = -1.0;
    
    return std::acos(cosTheta);
}

// For each particle, finds its nearest neighbor and returns:
// - The minimum angular distance to any other particle
// - The index of that nearest particle
std::vector<std::pair<double, size_t>> findNearestNeighbors(const std::vector<kdp::Vec3>& particles) {
    std::vector<std::pair<double, size_t>> nearestNeighbors;
    nearestNeighbors.reserve(particles.size());
    
    // For each particle, find its closest neighbor
    for (size_t i = 0; i < particles.size(); ++i) {
        double minDistance = M_PI;  // Start with maximum possible angular distance
        size_t nearestIdx = i;
        
        // Compare with all other particles
        for (size_t j = 0; j < particles.size(); ++j) {
            if (i != j) {
                double distance = angularDistance(particles[i], particles[j]);
                if (distance < minDistance) {
                    minDistance = distance;
                    nearestIdx = j;
                }
            }
        }
        
        nearestNeighbors.emplace_back(minDistance, nearestIdx);
    }
    
    return nearestNeighbors;
}

// Calculate the sum of squared energy fractions <f|f>
double calculateEnergyFractionSum(const std::vector<kdp::Vec3>& particles) {
    double ff = 0.0;
    auto particlesP = PowerSpectrum::PhatF::To_PhatF_Vec(particles);
    
    for(const auto& p : particlesP) {
        ff += kdp::Squared(p.f);
    }
    
    return ff;
}

// Reads particle data from a formatted input file
// Returns a vector of 4-vectors representing the particles
std::vector<kdp::Vec4> readParticlesFromFile(const std::string& filename) {
    std::vector<kdp::Vec4> particles;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open input file: " << filename << std::endl;
        return particles;
    }

    // Skip the header line
    std::string line;
    std::getline(file, line);

    // Process each line of particle data
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        int particle_num, id, status;
        double E, px, py, pz, mass, charge, visible;
        
        if (iss >> particle_num >> id >> status >> E >> px >> py >> pz >> mass >> charge >> visible) {
            // Create a 3-momentum vector
            kdp::Vec3 p3(px, py, pz);
            
            // Calculate proper relativistic energy: E = sqrt(p² + m²)
            double pMag = p3.Mag();
            double properE = std::sqrt(pMag*pMag + mass*mass);
            
            // Create 4-vector ensuring proper mass-shell condition
            particles.emplace_back(0., p3, kdp::Vec4from2::Mass);
        }
    }

    return particles;
}

int main(int argc, char* argv[]) {
    // Set up and verify configuration path
    const std::string configFilePath = "/home/mithila/EmPowerSpectrum/CONF/PowerSpectrum.conf";
    if (!std::ifstream(configFilePath).good()) {
        std::cerr << "Error: Configuration file " << configFilePath << " not found." << std::endl;
        return EXIT_FAILURE;
    }

    // Parse configuration using Qt settings
    QSettings parsedINI(configFilePath.c_str(), QSettings::NativeFormat);
    
    // Get required configuration parameters
    std::string dataFile = parsedINI.value("main/inputDataFile", "").toString().toStdString();
    if (dataFile.empty()) {
        std::cerr << "Error: No input data file specified in the configuration file." << std::endl;
        return EXIT_FAILURE;
    }
    
    size_t const lMax = size_t(parsedINI.value("power/lMax", 1024).toInt());
    std::string outputDir = parsedINI.value("main/outputDir", "./").toString().toStdString();
    if (outputDir.back() != '/' && outputDir.back() != '\\') {
        outputDir += "/";
    }

    // Get eta max values from config
    double etaMax_cal = parsedINI.value("detector/etaMax_cal", 3.0).toDouble();
    double etaMax_track = parsedINI.value("detector/etaMax_track", 3.0).toDouble();

    // Initialize detector and power spectrum computer
    ArrogantDetector_Lepton detector(parsedINI);
    PowerSpectrum Hl_computer(parsedINI.value("main/numThreads", 4).toInt());

    // Read particle data from input file
    std::vector<kdp::Vec4> inputParticles = readParticlesFromFile(dataFile);
    if (inputParticles.empty()) {
        std::cerr << "Error: No particles read from input file" << std::endl;
        return EXIT_FAILURE;
    }

    // Process particles through detector as both tracks and towers
    detector.Clear();
    detector(inputParticles, inputParticles);

    // Create container for power spectra (will hold 4 spectra)
    std::vector<std::vector<double>> Hl_set;

    // Variables to store energy fraction sums
    double tracks_ff = 0.0;
    double towers_ff = 0.0;

    // Process tracks and compute their power spectra
    {
        // Find nearest neighbors for Gaussian smoothing width
        auto track_NN = findNearestNeighbors(detector.Tracks());
        
        // Calculate mean nearest neighbor distance and energy fraction sum
        double mean = 0.;
        tracks_ff = calculateEnergyFractionSum(detector.Tracks());
        {
            auto tracks = PowerSpectrum::PhatF::To_PhatF_Vec(detector.Tracks());
            assert(track_NN.size() == tracks.size());
            
            for(size_t i = 0; i < track_NN.size(); ++i) {
                mean += track_NN[i].first * tracks[i].f;
            }
        }

        // Compute unsmoothed track spectrum
        PowerSpectrum::ShapedParticleContainer track_container(detector.Tracks());
        Hl_set.push_back(Hl_computer.Hl_Obs(lMax, track_container));

        // Compute smoothed track spectrum
        PowerSpectrum::ShapedParticleContainer track_container_smoothed(
            detector.Tracks(),
            h_Gaussian(std::sin(0.25*mean)*std::sqrt(-2./std::log1p(-0.9)))
        );
        Hl_set.push_back(Hl_computer.Hl_Obs(lMax, track_container_smoothed));
    }

    // Process towers and compute their power spectra
    {
        std::vector<kdp::Vec3> towers(detector.Towers().begin(), detector.Towers().end());
        auto tower_NN = findNearestNeighbors(towers);
        
        // Calculate mean nearest neighbor distance and energy fraction sum
        double mean = 0.;
        towers_ff = calculateEnergyFractionSum(towers);
        {
            auto towers_p = PowerSpectrum::PhatF::To_PhatF_Vec(towers);
            assert(tower_NN.size() == towers_p.size());
            
            for(size_t i = 0; i < tower_NN.size(); ++i) {
                mean += tower_NN[i].first * towers_p[i].f;
            }
        }

        // Compute unsmoothed tower spectrum
        PowerSpectrum::ShapedParticleContainer tower_container(towers);
        Hl_set.push_back(Hl_computer.Hl_Obs(lMax, tower_container));

        // Compute smoothed tower spectrum
        PowerSpectrum::ShapedParticleContainer tower_container_smoothed(
            towers,
            h_Gaussian(std::sin(0.25*mean)*std::sqrt(-2./std::log1p(-0.9)))
        );
        Hl_set.push_back(Hl_computer.Hl_Obs(lMax, tower_container_smoothed));
    }

    // Construct output filename for power spectrum data
    std::string inputFileName = GetProcess(dataFile);
    std::string outputFile = outputDir + "H_" + inputFileName + "_" + 
      std::to_string(int(std::round(kdp::ToDegrees(detector.GetSettings().squareWidth)))) + 
      "_" + std::to_string(etaMax_cal) + "_" + std::to_string(etaMax_track) + ".dat";
    
    // Write power spectrum data
    PowerSpectrum::WriteToFile(outputFile, Hl_set);
    
    // Create and write to log file with log_ prefix
    std::string logFile = outputDir + "log_" + inputFileName + "_" + 
      std::to_string(int(std::round(kdp::ToDegrees(detector.GetSettings().squareWidth)))) + 
      "_" + std::to_string(etaMax_cal) + "_" + std::to_string(etaMax_track) + ".dat";
    
    // Write summary to log file
    std::ofstream logStream(logFile);
    if (logStream.is_open()) {
      logStream << "Summary:" << std::endl
        << "Input file: " << dataFile << std::endl
        << "Output file: " << outputFile << std::endl
        << "Number of input particles: " << inputParticles.size() << std::endl
        << "Number of tracks: " << detector.Tracks().size() << std::endl
        << "Number of towers: " << detector.Towers().size() << std::endl
        << "Maximum l value: " << lMax << std::endl
        << "Square width: " << kdp::ToDegrees(detector.GetSettings().squareWidth) << " degrees" << std::endl
        << "Eta max (cal): " << etaMax_cal << std::endl
        << "Eta max (track): " << etaMax_track << std::endl
        << "Tracks <f|f>: " << tracks_ff << std::endl
        << "Towers <f|f>: " << towers_ff << std::endl;
      logStream.close();
    } else {
      std::cerr << "Error: Could not create log file: " << logFile << std::endl;
    }
    
    // Print summary of processing results
    std::cout << "\nProcessing summary:" << std::endl;
    std::cout << "Input file: " << dataFile << std::endl;
    std::cout << "Output file: " << outputFile << std::endl;
    std::cout << "Number of input particles: " << inputParticles.size() << std::endl;
    std::cout << "Number of tracks: " << detector.Tracks().size() << std::endl;
    std::cout << "Number of towers: " << detector.Towers().size() << std::endl;
    std::cout << "Maximum l value: " << lMax << std::endl;
    std::cout << "Square width: " << kdp::ToDegrees(detector.GetSettings().squareWidth) << " degrees" << std::endl;
    std::cout << "Eta max (cal): " << etaMax_cal << std::endl;
    std::cout << "Eta max (track): " << etaMax_track << std::endl;
    std::cout << "Tracks <f|f>: " << tracks_ff << std::endl;
    std::cout << "Towers <f|f>: " << towers_ff << std::endl;

    return 0;
}
