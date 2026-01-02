// Batch_H_computer.cpp
// Batch wrapper of H_computer.cpp: process every file in inputDir
// Behavior and outputs identical to H_computer.cpp, but loops over files in inputDir.

#include "Post_PYTHIA_PowerJets.hpp"
#include "PowerSpectrum.hpp"
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
#include <limits>
#include <QtCore/QSettings>
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>

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
        
        // Skip comment lines or empty lines
        if (line.find('#') != std::string::npos || line.empty()) {
            continue;
        }
        
        // Try to parse as first format first
        int particle_num, id, status;
        double E, px, py, pz, mass, charge, visible;
        
        // Read all fields including charge and visibility
        if (iss >> particle_num >> id >> status >> E >> px >> py >> pz >> mass >> charge >> visible) {
            kdp::Vector3<double> p3(px, py, pz);
            
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

// Helper: check if a path is a regular file
bool is_regular_file(const std::string &path) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0) return false;
    return S_ISREG(st.st_mode);
}

// Core per-file processing extracted from original H_computer.cpp main(). With added functionality as of 11-06-25.
// Keeps behavior unchanged; returns 0 on success, non-zero on failure for that file.
int ProcessSingleDataFile(const std::string& configFilePath, const std::string& dataFile) {
    // Validate data file
    if (!std::ifstream(dataFile).good()) {
        std::cerr << "Error: Input data file does not exist: " << dataFile << std::endl;
        return EXIT_FAILURE;
    }

    // QSettings read
    QSettings parsedINI(configFilePath.c_str(), QSettings::NativeFormat);

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
    
    std::cout << "Configuration for file: " << dataFile << std::endl;
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
    
    // Create modified configuration for detector (to change etaMax and beam_holes if needed)
    std::string modifiedConfigPath = configFilePath + ".tmp";
    
    {
        std::ifstream originalConfig(configFilePath);
        std::ofstream modifiedConfig(modifiedConfigPath);
        
        if (originalConfig.is_open() && modifiedConfig.is_open()) {
            std::string line;
            while (std::getline(originalConfig, line)) {
                if (line.find("etaMax_cal") != std::string::npos) {
                    modifiedConfig << "etaMax_cal = " << effectiveEtaMax_cal << std::endl;
                } else if (line.find("etaMax_track") != std::string::npos) {
                    modifiedConfig << "etaMax_track = " << effectiveEtaMax_track << std::endl;
                } else if (line.find("beam_holes") != std::string::npos) {
                    modifiedConfig << "beam_holes = " << (useBeamHoles ? "ON" : "OFF") << std::endl;
                } else {
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

	// NEW: make sure the per-file temporary config contains the input filename
	{
	  QSettings tmpSetter(modifiedConfigPath.c_str(), QSettings::NativeFormat);
	  tmpSetter.setValue("main/inputDataFile", QString::fromStdString(dataFile));
	  tmpSetter.sync();
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
      for (const auto& particle : particleData) {
        if (particle.isVisible) {
          filteredParticles.push_back(particle.momentum);
        }
      }
      std::cout << "Using all " << filteredParticles.size() << " visible particles (beam_holes=OFF)" << std::endl;
    }

    // Set filtered particles for power spectrum calculation
    powerJets.SetParticles(filteredParticles);    
    
    // Extract the file name without the extension using the existing GetProcess helper
    std::string inputFileName = GetProcess(dataFile);

    // Retrieve output directory from the configuration file
    std::string outputDir = parsedINI.value("main/outputDir", "./").toString().toStdString();

    // Ensure output directory ends with a slash
    if (!outputDir.empty() && outputDir.back() != '/' && outputDir.back() != '\\') {
        outputDir += "/";
    }

    // Check if the output directory exists; if not, use the current directory
    if (!outputDir.empty() && !std::ifstream(outputDir).good()) {
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

    // Extra: store detector-level directions & square width (for AngularResolution)
    std::vector<kdp::Vector3<double>> detectorTracks;
    std::vector<kdp::Vector3<double>> detectorTowers;   // will exclude the missing-E tower
    double squareWidthRadians = 0.0;


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

            // Store detector-level tracks/towers (exclude missing-E tower) and square width in radians
            detectorTracks.assign(detector->Tracks().begin(), detector->Tracks().end());
            if (!detector->Towers().empty()) {
                detectorTowers.assign(detector->Towers().begin(), detector->Towers().end() - 1);
            }
            squareWidthRadians = detector->GetSettings().squareWidth;

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

            // Store detector-level tracks/towers (exclude missing-E tower) and square width in radians
            detectorTracks.assign(completeDetector->Tracks().begin(),
                                  completeDetector->Tracks().end());
            if (!completeDetector->Towers().empty()) {
                detectorTowers.assign(completeDetector->Towers().begin(),
                                      completeDetector->Towers().end() - 1);
            }
            squareWidthRadians = completeDetector->GetSettings().squareWidth;

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

    
    // ------------------------------------------------------------
    // Compute three angular scales using only objects within |eta| <= etaCut,
    // where etaCut = min(etaMax_cal, etaMax_track) from the CONF file.
    // This does NOT affect Hl computation; it's only for the angular scales.
    // ------------------------------------------------------------
    
    {
      using real_t = PowerSpectrum::real_t;
      
      // 0) Read eta limits from the original configuration
      const double etaMax_cal   = parsedINI.value("detector/etaMax_cal",   4.9).toDouble();
      const double etaMax_track = parsedINI.value("detector/etaMax_track", 4.9).toDouble();
      const double etaCut       = std::min(etaMax_cal, etaMax_track);
      
      // 1) Build the eta-cut subsets
      
      // 1a) Final-state (delta) set: from filteredParticles (same set used for Hl_FinalState)
      std::vector<const kdp::Vector3<double>*> particleDirs_etaCut;
      particleDirs_etaCut.reserve(filteredParticles.size());
      for (const auto& fp : filteredParticles) {
	const kdp::Vector3<double>& p3 = fp.p(); // from Vector4::p()
	const double eta = p3.Eta();
	if (std::isfinite(eta) && std::fabs(eta) <= etaCut) {
	  particleDirs_etaCut.push_back(&p3);
	}
      }
      
      // 1b) Detector tracks/towers subsets (exclude missing-E tower as you already do)
      std::vector<kdp::Vector3<double>> tracks_etaCut;
      std::vector<kdp::Vector3<double>> towers_etaCut;
      
      tracks_etaCut.reserve(detectorTracks.size());
      for (const auto& v : detectorTracks) {
	const double eta = v.Eta();
	if (std::isfinite(eta) && std::fabs(eta) <= etaCut) {
	  tracks_etaCut.push_back(v);
	}
      }
      
      towers_etaCut.reserve(detectorTowers.size());
      for (const auto& v : detectorTowers) {
	const double eta = v.Eta();
	if (std::isfinite(eta) && std::fabs(eta) <= etaCut) {
	  towers_etaCut.push_back(v);
	}
      }
      
      // 2) ksi_min_delta: smallest angular separation within the eta-cut final-state set
      double ksi_min_delta = std::numeric_limits<double>::infinity();
      if (particleDirs_etaCut.size() >= 2) {
	for (size_t i = 0; i + 1 < particleDirs_etaCut.size(); ++i) {
	  const kdp::Vector3<double>& pi = *particleDirs_etaCut[i];
	  for (size_t j = i + 1; j < particleDirs_etaCut.size(); ++j) {
	    const kdp::Vector3<double>& pj = *particleDirs_etaCut[j];
	    const double angle = pi.InteriorAngle(pj);
	    if (angle < ksi_min_delta) {
	      ksi_min_delta = angle;
	    }
	  }
	}
      }
      
      // 3) Build PhatF containers for AngularResolution using a COMMON normalization
      //    to the subset's detected total energy (so Σ f over tracks+towers == 1).
      double subsetDetectedEnergy = 0.0;
      for (const auto& v : tracks_etaCut) {
	subsetDetectedEnergy += v.Mag();
      }
      for (const auto& v : towers_etaCut) {
	subsetDetectedEnergy += v.Mag();
      }
      
      std::vector<PowerSpectrum::PhatF> trackPhats_etaCut;
      std::vector<PowerSpectrum::PhatF> towerPhats_etaCut;
      if (subsetDetectedEnergy > 0.0) {
	const real_t totalE_in = static_cast<real_t>(subsetDetectedEnergy);
	trackPhats_etaCut = PowerSpectrum::PhatF::To_PhatF_Vec(tracks_etaCut, totalE_in);
	towerPhats_etaCut = PowerSpectrum::PhatF::To_PhatF_Vec(towers_etaCut, totalE_in);
      }
      
      // 4) Angular resolutions (as-is, using PowerSpectrum’s own logic)
      
      // 4a) Track-only effective resolution
      real_t ksi_track_only = PowerSpectrum::AngularResolution(trackPhats_etaCut);
      
      // 4b) Full sample resolution (tracks + towers + squareWidth)
      const real_t squareWidth_r = static_cast<real_t>(squareWidthRadians);
      real_t ksi_sample = PowerSpectrum::AngularResolution(
							   trackPhats_etaCut,
							   towerPhats_etaCut,
							   squareWidth_r
							   );
      
      // 5) Convert to lMax values: lMax = 2*pi / ksi (guard zero/inf) and round
      const double twoPi = 2.0 * M_PI;
      
      double lMax_delta      = (std::isfinite(ksi_min_delta) && ksi_min_delta > 0.0)
	? (twoPi / ksi_min_delta)
	: 0.0;
      
      double lMax_track_only = (std::isfinite(ksi_track_only) && ksi_track_only > real_t(0))
	? (twoPi / static_cast<double>(ksi_track_only))
	: 0.0;
      
      double lMax_smeared    = (std::isfinite(ksi_sample) && ksi_sample > real_t(0))
	? (twoPi / static_cast<double>(ksi_sample))
	: 0.0;
      
      const long lMax_delta_int      = static_cast<long>(std::round(lMax_delta));
      const long lMax_track_only_int = static_cast<long>(std::round(lMax_track_only));
      const long lMax_smeared_int    = static_cast<long>(std::round(lMax_smeared));
      
      // 6) Prepend the new header line to the existing .dat file (leave everything else intact)
      std::ifstream inFile(outputFilePath.c_str());
      if (!inFile.is_open()) {
	std::cerr << "Warning: Unable to reopen output file for header update: "
		  << outputFilePath << std::endl;
      } else {
	std::stringstream contents;
	contents << inFile.rdbuf();
	inFile.close();
	
	std::ofstream outFile(outputFilePath.c_str());
	if (!outFile.is_open()) {
	  std::cerr << "Warning: Unable to write updated output file: "
		    << outputFilePath << std::endl;
	} else {
	  outFile << "# lMax [delta, track-only-eff, sample-res]: ["
		  << lMax_delta_int << ", "
		  << lMax_track_only_int << ", "
		  << lMax_smeared_int << "]\n"
		  << contents.str();
	}
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
        
        logFile << "Square width: " << squareWidthDegrees << " degrees" << std::endl;
        logFile << "Detector type: " << detectorType.toStdString() << std::endl;
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
        
        double detectionEfficiency = (detectedTotalEnergy / inputTotalEnergy) * 100.0;
        logFile << "Detection efficiency: " << detectionEfficiency << "%" << std::endl;
        
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
    std::cout << "\nProcessing summary for file:" << std::endl;
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
    
    // Clean up temporary configuration file
    if (modifiedConfigPath != configFilePath) {
        std::remove(modifiedConfigPath.c_str());
    }

    return 0;
}

int main(int /*argc*/, char* /*argv*/[]) {
    // Set the path to the configuration file
    const std::string configFilePath = "/home/mithila/EmPowerSpectrum/CONF/Batch_PowerSpectrum.conf";

    // Check if the configuration file exists
    if (!std::ifstream(configFilePath).good()) {
        std::cerr << "Error: Configuration file " << configFilePath << " not found." << std::endl;
        return EXIT_FAILURE;
    }

    // Read config once here to obtain inputDir
    QSettings parsedINI(configFilePath.c_str(), QSettings::NativeFormat);

    // Prefer inputDir for batch mode. If not set, behave like original and try inputDataFile.
    std::string inputDir = parsedINI.value("main/inputDir", "").toString().toStdString();
    std::string singleInputFile = parsedINI.value("main/inputDataFile", "").toString().toStdString();

    if (inputDir.empty() && singleInputFile.empty()) {
        std::cerr << "Error: Neither main/inputDir nor main/inputDataFile is set in the configuration file." << std::endl;
        return EXIT_FAILURE;
    }

    // If inputDir is empty, fall back to single file (run it once)
    if (inputDir.empty()) {
        std::cout << "No inputDir found; processing single file from inputDataFile: " << singleInputFile << std::endl;
        return ProcessSingleDataFile(configFilePath, singleInputFile);
    }

    // Ensure inputDir ends with '/'
    if (inputDir.back() != '/' && inputDir.back() != '\\') inputDir += "/";

    // Try to open directory
    DIR* dir = opendir(inputDir.c_str());
    if (!dir) {
        std::cerr << "Error: Could not open input directory '" << inputDir << "': " << strerror(errno) << std::endl;
        return EXIT_FAILURE;
    }

    struct dirent* entry;
    int failCount = 0;
    int successCount = 0;
    while ((entry = readdir(dir)) != nullptr) {
        std::string name(entry->d_name);
        if (name == "." || name == "..") continue;
        if (name.size() > 0 && name[0] == '.') continue; // skip hidden files

        std::string fullpath = inputDir + name;

        // Only process regular files
        if (!is_regular_file(fullpath)) {
            std::cerr << "Skipping non-regular file: " << fullpath << std::endl;
            continue;
        }

        std::cout << "\n=== Processing file: " << fullpath << " ===\n";
        int rc = ProcessSingleDataFile(configFilePath, fullpath);
        if (rc == 0) {
            ++successCount;
        } else {
            ++failCount;
            std::cerr << "Processing failed for file: " << fullpath << " (rc=" << rc << ")\n";
        }
    }

    closedir(dir);

    std::cout << "\nBatch processing complete. Success: " << successCount << ", Failures: " << failCount << std::endl;

    return (failCount == 0) ? 0 : EXIT_FAILURE;
}
