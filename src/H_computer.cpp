#include "Post_PYTHIA_PowerJets.hpp"
#include "Admin.hpp"
#include "ArrogantDetector.hpp"
#include "CompleteArrogantDetector.hpp"
#include "kdpVectors.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <QtCore/QSettings>

// Structure to hold particle data with its charge information
struct ParticleData {
    kdp::Vector4<double> momentum;
    bool isCharged;
    bool isVisible;
};

// Reads particle data from a formatted input file
std::vector<ParticleData> readParticlesFromFile(const std::string& filename) {
    std::vector<ParticleData> particles;
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
        
        // Check if this is the first format (with Particle, ID, Status columns)
        // or the second format (with # pdgID columns)
        if (line.find('#') != std::string::npos || line.empty()) {
            // Skip comment lines or empty lines
            continue;
        }
        
        // Try to parse as first format first
        int particle_num, id, status;
        double E, px, py, pz, mass, charge, visible;
        
        // Read all fields including charge and visibility
        if (iss >> particle_num >> id >> status >> E >> px >> py >> pz >> mass >> charge >> visible) {
            // First format: Particle ID Status E px py pz m Charge Visible
            // Create a 3-momentum vector
            kdp::Vector3<double> p3(px, py, pz);
            
            // Create 4-vector with proper energy
            ParticleData particle;
            particle.momentum = kdp::Vector4<double>(E, p3, kdp::Vec4from2::Energy);
            particle.isCharged = (charge != 0);
            particle.isVisible = (visible > 0.5);
            
            particles.push_back(particle);
        } else {
            // Try second format: pdgID pt eta phi mass px py pz E Charge
            iss.clear();
            iss.str(line);
            
            int pdgID;
            double pt, eta, phi, mass2, px2, py2, pz2, E2, charge2;
            
            if (iss >> pdgID >> pt >> eta >> phi >> mass2 >> px2 >> py2 >> pz2 >> E2 >> charge2) {
                kdp::Vector3<double> p3(px2, py2, pz2);
                
                ParticleData particle;
                // Use the constructor that takes energy and 3-vector
                particle.momentum = kdp::Vector4<double>(E2, p3, kdp::Vec4from2::Energy);
                particle.isCharged = (charge2 != 0);
                particle.isVisible = true;  // Assume visible in second format
                
                particles.push_back(particle);
            }
        }
    }

    return particles;
}

// Calculate total energy of particles
double calculateTotalEnergy(const std::vector<ParticleData>& particles) {
    double totalEnergy = 0.0;
    for (const auto& particle : particles) {
        totalEnergy += particle.momentum.x0;  // x0 is the energy component
    }
    return totalEnergy;
}

// Custom function to analyze particle distribution based on eta values
void analyzeParticleDistribution(
    const std::vector<ParticleData>& particles, 
    double etaMax,
    double& energyWithinEta,
    double& energyOutsideEta,
    int& countWithinEta,
    int& countOutsideEta) {
    
    energyWithinEta = 0.0;
    energyOutsideEta = 0.0;
    countWithinEta = 0;
    countOutsideEta = 0;
    
    for (const auto& particle : particles) {
        double eta = std::fabs(particle.momentum.p().Eta());
        double energy = particle.momentum.x0;
        
        if (eta <= etaMax) {
            energyWithinEta += energy;
            countWithinEta++;
        } else {
            energyOutsideEta += energy;
            countOutsideEta++;
        }
    }
}

int main(int /*argc*/, char* /*argv*/[]) {
    // Set the path to the configuration file
    const std::string configFilePath = "/home/mithila/EmPowerSpectrum/CONF/PowerSpectrum.conf";

    // Check if the configuration file exists
    if (!std::ifstream(configFilePath).good()) {
        std::cerr << "Error: Configuration file " << configFilePath << " not found." << std::endl;
        return EXIT_FAILURE;
    }

    // Setup QSettings with the correct path to the configuration file
    QSettings parsedINI(configFilePath.c_str(), QSettings::NativeFormat);

    // Retrieve the input data file from the configuration file
    std::string dataFile = parsedINI.value("main/inputDataFile", "").toString().toStdString();
    if (dataFile.empty()) {
        std::cerr << "Error: No input data file specified in the configuration file." << std::endl;
        return EXIT_FAILURE;
    }
    
    // Check if the input file exists before proceeding
    if (!std::ifstream(dataFile).good()) {
        std::cerr << "Error: Input data file does not exist: " << dataFile << std::endl;
        return EXIT_FAILURE;
    }

    // Read beam_holes setting (ON/OFF)
    QString beamHolesSetting = parsedINI.value("detector/beam_holes", "ON").toString().toUpper();
    bool useBeamHoles = (beamHolesSetting == "ON" || beamHolesSetting == "TRUE" || beamHolesSetting == "1");
    
    // Store original eta values
    double originalEtaMax_cal = parsedINI.value("detector/etaMax_cal", 4.9).toDouble();
    double originalEtaMax_track = parsedINI.value("detector/etaMax_track", 4.9).toDouble();
    
    // Define extremely high etaMax value to use when beam_holes is OFF (for full spherical coverage)
    const double HIGH_ETA_MAX = 999.0;  // Effectively infinite
    
    // If beam_holes is OFF, use the high eta max for processing
    double effectiveEtaMax_cal = useBeamHoles ? originalEtaMax_cal : HIGH_ETA_MAX;
    double effectiveEtaMax_track = useBeamHoles ? originalEtaMax_track : HIGH_ETA_MAX;
    
    // Get detector type from configuration
    QString detectorType = parsedINI.value("detector/type", "lepton").toString().toLower();
    
    std::cout << "Configuration:" << std::endl;
    std::cout << "  Beam holes: " << (useBeamHoles ? "ON" : "OFF") << std::endl;
    std::cout << "  Detector type: " << detectorType.toStdString() << std::endl;
    std::cout << "  Original Eta max (cal): " << originalEtaMax_cal << std::endl;
    std::cout << "  Original Eta max (track): " << originalEtaMax_track << std::endl;
    if (!useBeamHoles) {
        std::cout << "  Effective Eta max: " << HIGH_ETA_MAX << " (beam_holes=OFF)" << std::endl;
    }
    
    // Read the raw particle data with charge information
    std::vector<ParticleData> particleData = readParticlesFromFile(dataFile);
    
    if (particleData.empty()) {
        std::cerr << "Error: No particles read from input file: " << dataFile << std::endl;
        return EXIT_FAILURE;
    }
    
    // Create modified configuration for detector
    std::string modifiedConfigPath = configFilePath + ".tmp";
    
    {
        std::ifstream originalConfig(configFilePath);
        std::ofstream modifiedConfig(modifiedConfigPath);
        
        if (originalConfig.is_open() && modifiedConfig.is_open()) {
            std::string line;
            while (std::getline(originalConfig, line)) {
                // Find and replace etaMax values
                if (line.find("etaMax_cal") != std::string::npos) {
                    modifiedConfig << "etaMax_cal = " << effectiveEtaMax_cal << std::endl;
                } else if (line.find("etaMax_track") != std::string::npos) {
                    modifiedConfig << "etaMax_track = " << effectiveEtaMax_track << std::endl;
                } else if (line.find("beam_holes") != std::string::npos) {
                    // Replace with our actual setting
                    modifiedConfig << "beam_holes = " << (useBeamHoles ? "ON" : "OFF") << std::endl;
                } else {
                    // Keep other lines unchanged, including detector type
                    modifiedConfig << line << std::endl;
                }
            }
            
            originalConfig.close();
            modifiedConfig.close();
            
            if (!useBeamHoles) {
                std::cout << "Created modified configuration with etaMax = " << HIGH_ETA_MAX 
                        << " to capture all particles" << std::endl;
            } else {
                std::cout << "Created modified configuration with original etaMax values" << std::endl;
            }
        } else {
            std::cerr << "Warning: Could not create modified configuration. Using original values." << std::endl;
            modifiedConfigPath = configFilePath;  // Fall back to original config
        }
    }

    // Separate particles into charged, neutral, and invisible vectors
    std::vector<kdp::Vector4<double>> chargedParticles;
    std::vector<kdp::Vector4<double>> neutralParticles;
    std::vector<kdp::Vector4<double>> invisibleParticles;
    
    for (const auto& particle : particleData) {
        if (!particle.isVisible) {
            invisibleParticles.push_back(particle.momentum);
        } else if (particle.isCharged) {
            chargedParticles.push_back(particle.momentum);
        } else {
            neutralParticles.push_back(particle.momentum);
        }
    }
    
    // Create a vector of just momenta for total energy calculation
    std::vector<kdp::Vector4<double>> allMomenta;
    for (const auto& particle : particleData) {
        allMomenta.push_back(particle.momentum);
    }

    // Initialize Post_PYTHIA_PowerJets with the modified configuration
    Post_PYTHIA_PowerJets powerJets(modifiedConfigPath);
    if (powerJets.GetStatus() != Post_PYTHIA_PowerJets::Status::OK) {
        std::cerr << "Error: Failed to initialize Post_PYTHIA_PowerJets properly." << std::endl;
        if (modifiedConfigPath != configFilePath) {
            std::remove(modifiedConfigPath.c_str());  // Clean up temp file
        }
        return EXIT_FAILURE;
    }

    // Filter particles based on beam_holes setting and pass to Post_PYTHIA_PowerJets
    std::vector<kdp::Vector4<double>> filteredParticles;
    if (useBeamHoles) {
      // Filter particles within eta range
      for (const auto& particle : particleData) {
        if (particle.isVisible) {
	  double eta = std::fabs(particle.momentum.p().Eta());
	  if (eta <= originalEtaMax_cal) {
	    filteredParticles.push_back(particle.momentum);
	  }
        }
      }
      std::cout << "Filtered to " << filteredParticles.size() << " visible particles within eta <= " << originalEtaMax_cal << std::endl;
    } else {
      // Use all visible particles when beam_holes is OFF
      for (const auto& particle : particleData) {
        if (particle.isVisible) {
	  filteredParticles.push_back(particle.momentum);
        }
      }
      std::cout << "Using all " << filteredParticles.size() << " visible particles (beam_holes=OFF)" << std::endl;
    }

    // Set filtered particles for power spectrum calculation
    powerJets.SetParticles(filteredParticles);    
    
    // Extract the file name without the extension
    std::string inputFileName = GetProcess(dataFile);

    // Retrieve output directory from the configuration file
    std::string outputDir = parsedINI.value("main/outputDir", "./").toString().toStdString();

    // Ensure output directory ends with a slash
    if (outputDir.back() != '/' && outputDir.back() != '\\') {
        outputDir += "/";
    }

    // Check if the output directory exists; if not, use the current directory
    if (!std::ifstream(outputDir).good()) {
        std::cerr << "Warning: Output directory does not exist: " << outputDir << std::endl;
        outputDir = "./";
    }

    // Construct the output file name (use same improved naming as the log file)
    std::string outputFilePath;
    if (useBeamHoles) {
      outputFilePath = outputDir + "H_" + inputFileName + "_with_holes_" +
	std::to_string(originalEtaMax_cal) + ".dat";
    } else {
      outputFilePath = outputDir + "H_" + inputFileName + "_spherical.dat";
    }
 
    // Write the power spectrum to the output file using Post_PYTHIA_PowerJets
    powerJets.WritePowerSpectrum(outputFilePath);
    
    // Calculate total energy of input particles
    double inputTotalEnergy = calculateTotalEnergy(particleData);
    
    // Analyze particle distribution before detector processing
    double energyWithinEta, energyOutsideEta;
    int countWithinEta, countOutsideEta;
    analyzeParticleDistribution(
        particleData, originalEtaMax_cal, 
        energyWithinEta, energyOutsideEta, 
        countWithinEta, countOutsideEta);
    
    // Initialize detector settings
    QSettings modifiedSettings(modifiedConfigPath.c_str(), QSettings::NativeFormat);
    
    // Variables to store detector output
    int numTracks = 0;
    int numTowers = 0;
    double trackEnergy = 0.0;
    double towerEnergy = 0.0;
    double missingEnergy = 0.0;
    double detectedTotalEnergy = 0.0;
    double squareWidthDegrees = 0.0;  // Store square width directly instead of settings object

    // Initialize the appropriate detector based on configuration
    if (useBeamHoles) {
        // Use regular ArrogantDetector with beam holes when beam_holes=ON
        ArrogantDetector* detector = nullptr;
        try {
            if (detectorType == "hadron") {
                detector = new ArrogantDetector_Hadron(modifiedSettings);
            } else {
                detector = new ArrogantDetector_Lepton(modifiedSettings);
            }
            
            // Process particles through the detector with beam holes
            detector->Clear();
            (*detector)(neutralParticles, chargedParticles, invisibleParticles);
            
            // Get detector data
            numTracks = static_cast<int>(detector->Tracks().size());
            numTowers = static_cast<int>(detector->Towers().size());
            
            // Calculate energy in tracks
            for (const auto& track : detector->Tracks()) {
                double energy = track.Mag();
                trackEnergy += energy;
                detectedTotalEnergy += energy;
            }
            
            // Calculate energy in towers (excluding missing energy tower)
            for (size_t i = 0; i < detector->Towers().size() - 1; i++) {
                const auto& tower = detector->Towers()[i];
                double energy = tower.Mag();
                towerEnergy += energy;
                detectedTotalEnergy += energy;
            }
            
            // Get missing energy (the last tower)
            if (!detector->Towers().empty()) {
                missingEnergy = detector->Towers().back().Mag();
            }
            
            // Get square width in degrees directly
            squareWidthDegrees = kdp::ToDegrees(detector->GetSettings().squareWidth);

            // Clean up
            delete detector;
        }
        catch (const std::exception& e) {
            std::cerr << "Error initializing detector: " << e.what() << std::endl;
            if (modifiedConfigPath != configFilePath) {
                std::remove(modifiedConfigPath.c_str());  // Clean up temp file
            }
            return EXIT_FAILURE;
        }
    } 
    else {
        // Use CompleteArrogantDetector with no beam holes when beam_holes=OFF
        CompleteArrogantDetector* completeDetector = nullptr;
        try {
            if (detectorType == "hadron") {
                completeDetector = new CompleteArrogantDetector_Hadron(modifiedSettings);
            } else {
                completeDetector = new CompleteArrogantDetector_Lepton(modifiedSettings);
            }
            
            // Process particles through the complete detector without beam holes
            completeDetector->Clear();
            
            // Now call the detector with properly separated particles
            (*completeDetector)(neutralParticles, chargedParticles, invisibleParticles);
            
            // Get detector data
            numTracks = static_cast<int>(completeDetector->Tracks().size());
            numTowers = static_cast<int>(completeDetector->Towers().size());
            
            // Calculate energy in tracks
            for (const auto& track : completeDetector->Tracks()) {
                double energy = track.Mag();
                trackEnergy += energy;
                detectedTotalEnergy += energy;
            }
            
            // Calculate energy in towers (excluding missing energy tower)
            for (size_t i = 0; i < completeDetector->Towers().size() - 1; i++) {
                const auto& tower = completeDetector->Towers()[i];
                double energy = tower.Mag();
                towerEnergy += energy;
                detectedTotalEnergy += energy;
            }
            
            // Get missing energy (the last tower)
            if (!completeDetector->Towers().empty()) {
                missingEnergy = completeDetector->Towers().back().Mag();
            }
            
            // Get square width in degrees directly
            squareWidthDegrees = kdp::ToDegrees(completeDetector->GetSettings().squareWidth);

            // Clean up
            delete completeDetector;
        }
        catch (const std::exception& e) {
            std::cerr << "Error initializing complete detector: " << e.what() << std::endl;
            if (modifiedConfigPath != configFilePath) {
                std::remove(modifiedConfigPath.c_str());  // Clean up temp file
            }
            return EXIT_FAILURE;
        }
    }

    // Get lMax from configuration
    int lMax = parsedINI.value("power/lMax", 256).toInt();
    
    // Create the log file with improved naming
    std::string logFilePath;
    if (useBeamHoles) {
        logFilePath = outputDir + "log_" + inputFileName + "_with_holes_" + 
                   std::to_string(originalEtaMax_cal) + ".txt";
    } else {
        logFilePath = outputDir + "log_" + inputFileName + "_spherical.txt";
    }
                           
    std::ofstream logFile(logFilePath);
    
    if (logFile.is_open()) {
        // Write the information to the log file
        logFile << "Number of tracks: " << numTracks << std::endl;
        logFile << "Number of towers: " << numTowers << std::endl;
        logFile << "Maximum l value: " << lMax << std::endl;
        
        // Write square width in degrees
        logFile << "Square width: " << squareWidthDegrees << " degrees" << std::endl;
        
        // Write detector type
        logFile << "Detector type: " << detectorType.toStdString() << std::endl;
        
        // Write eta max values and beam holes settings
        logFile << "Original Eta max (cal): " << originalEtaMax_cal << std::endl;
        logFile << "Original Eta max (track): " << originalEtaMax_track << std::endl;
        if (!useBeamHoles) {
            logFile << "Effective Eta max: " << HIGH_ETA_MAX << " (beam_holes=OFF)" << std::endl;
            logFile << "Using CompleteArrogantDetector (truly spherical detector)" << std::endl;
        } else {
            logFile << "Using ArrogantDetector (with beam holes)" << std::endl;
        }
        logFile << "Beam holes: " << (useBeamHoles ? "ON" : "OFF") << std::endl;
        
        // Add particle distribution info
        logFile << "Particles within eta " << originalEtaMax_cal << ": " << countWithinEta << std::endl;
        logFile << "Particles outside eta " << originalEtaMax_cal << ": " << countOutsideEta << std::endl;
        logFile << "Energy within eta " << originalEtaMax_cal << ": " << energyWithinEta 
                << " GeV (" << (energyWithinEta * 100.0 / inputTotalEnergy) << "%)" << std::endl;
        logFile << "Energy outside eta " << originalEtaMax_cal << ": " << energyOutsideEta 
                << " GeV (" << (energyOutsideEta * 100.0 / inputTotalEnergy) << "%)" << std::endl;
        
        // Add energy information
        logFile << "Total energy of input particles: " << inputTotalEnergy << " GeV" << std::endl;
        logFile << "Energy in tracks: " << trackEnergy << " GeV" << std::endl;
        logFile << "Energy in calorimeter towers: " << towerEnergy << " GeV" << std::endl;
        logFile << "Missing energy (last tower): " << missingEnergy << " GeV" << std::endl;
        logFile << "Total energy of detected particles: " << detectedTotalEnergy << " GeV" << std::endl;
        
        // Calculate and record detection efficiency
        double detectionEfficiency = (detectedTotalEnergy / inputTotalEnergy) * 100.0;
        logFile << "Detection efficiency: " << detectionEfficiency << "%" << std::endl;
        
        // Calculate and record missing energy percentage
        double missingEnergyPercentage = (missingEnergy / inputTotalEnergy) * 100.0;
        logFile << "Missing energy percentage: " << missingEnergyPercentage << "%" << std::endl;
        
        // Additional information for diagnosing detector coverage
        logFile << "\nDetector Coverage Analysis:" << std::endl;
        logFile << "  Input particle count: " << particleData.size() << std::endl;
        logFile << "  Charged particles: " << chargedParticles.size() << std::endl;
        logFile << "  Neutral particles: " << neutralParticles.size() << std::endl;
        logFile << "  Invisible particles: " << invisibleParticles.size() << std::endl;
        
        logFile.close();
        std::cout << "Log file created at: " << logFilePath << std::endl;
    } else {
        std::cerr << "Warning: Unable to create log file at " << logFilePath << std::endl;
    }
    
    // Print summary to console
    std::cout << "\nProcessing summary:" << std::endl;
    std::cout << "Input file: " << dataFile << std::endl;
    std::cout << "Output file: " << outputFilePath << std::endl;
    std::cout << "Number of tracks: " << numTracks << std::endl;
    std::cout << "Number of towers: " << numTowers << std::endl;
    std::cout << "Maximum l value: " << lMax << std::endl;
    std::cout << "Square width: " << squareWidthDegrees << " degrees" << std::endl;
    std::cout << "Detector type: " << detectorType.toStdString() << std::endl;
    if (!useBeamHoles) {
        std::cout << "Using: CompleteArrogantDetector (truly spherical detector)" << std::endl;
    } else {
        std::cout << "Using: ArrogantDetector (with beam holes)" << std::endl;
    }
    std::cout << "Original Eta max (cal): " << originalEtaMax_cal << std::endl;
    std::cout << "Original Eta max (track): " << originalEtaMax_track << std::endl;
    if (!useBeamHoles) {
        std::cout << "Effective Eta max: " << HIGH_ETA_MAX << " (beam_holes=OFF)" << std::endl;
    }
    std::cout << "Beam holes: " << (useBeamHoles ? "ON" : "OFF") << std::endl;
    
    // Add particle counts
    std::cout << "Particles within eta " << originalEtaMax_cal << ": " << countWithinEta << std::endl;
    std::cout << "Particles outside eta " << originalEtaMax_cal << ": " << countOutsideEta << std::endl;
    std::cout << "Energy within eta " << originalEtaMax_cal << ": " << energyWithinEta 
              << " GeV (" << (energyWithinEta * 100.0 / inputTotalEnergy) << "%)" << std::endl;
    std::cout << "Energy outside eta " << originalEtaMax_cal << ": " << energyOutsideEta 
              << " GeV (" << (energyOutsideEta * 100.0 / inputTotalEnergy) << "%)" << std::endl;
    
    std::cout << "Total energy of input particles: " << inputTotalEnergy << " GeV" << std::endl;
    std::cout << "Energy in tracks: " << trackEnergy << " GeV" << std::endl;
    std::cout << "Energy in calorimeter towers: " << towerEnergy << " GeV" << std::endl;
    std::cout << "Missing energy (last tower): " << missingEnergy << " GeV" << std::endl;
    std::cout << "Total energy of detected particles: " << detectedTotalEnergy << " GeV" << std::endl;
    
    // Calculate and display detection efficiency
    double detectionEfficiency = (detectedTotalEnergy / inputTotalEnergy) * 100.0;
    std::cout << "Detection efficiency: " << detectionEfficiency << "%" << std::endl;
    
    // Calculate and display missing energy percentage
    double missingEnergyPercentage = (missingEnergy / inputTotalEnergy) * 100.0;
    std::cout << "Missing energy percentage: " << missingEnergyPercentage << "%" << std::endl;
    
    // Print particle type counts
    std::cout << "\nParticle type counts:" << std::endl;
    std::cout << "  Charged particles: " << chargedParticles.size() << std::endl;
    std::cout << "  Neutral particles: " << neutralParticles.size() << std::endl;
    std::cout << "  Invisible particles: " << invisibleParticles.size() << std::endl;
    
    // Clean up temporary configuration file
    if (modifiedConfigPath != configFilePath) {
        std::remove(modifiedConfigPath.c_str());
    }

    return 0;
}
