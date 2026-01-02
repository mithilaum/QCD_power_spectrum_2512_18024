#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <iomanip>
#include <string>
#include <sys/stat.h>
#include "VectorGenerators.hpp"
#include "ConfigParser.hpp"

// Function to create directory if it doesn't exist
bool createDirectory(const std::string& path) {
#ifdef _WIN32
    return _mkdir(path.c_str()) == 0 || errno == EEXIST;
#else
    return mkdir(path.c_str(), 0777) == 0 || errno == EEXIST;
#endif
}

// Function to calculate momentum magnitude from fixed energy and mass
double calculateMomentumMagnitude(double energy, double mass) {
    // Check if energy is sufficient for mass-shell condition
    double energySquared = energy * energy;
    double massSquared = mass * mass;
    
    if (energySquared <= massSquared) {
        throw std::runtime_error("Mass-shell condition cannot be satisfied: Energy (" + 
                               std::to_string(energy) + " GeV) is not greater than mass (" + 
                               std::to_string(mass) + " GeV)");
    }
    
    return std::sqrt(energySquared - massSquared);
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file>" << std::endl;
        return 1;
    }

    try {
        // Parse configuration
        ConfigParser config(argv[1]);
        
        // Get configuration values
        int N = config.getInt("num_vectors");
        std::string method = config.getString("generation_method", "fibonacci");
        double m = config.getDouble("mass", 1.3957000000000000e-01);
        double energy = config.getDouble("energy", 1.0);  // Default energy of 1.0 GeV
        int status = config.getInt("status", 83);
        int ID = config.getInt("particle_id", 211);
        int base_charge = config.getInt("charge", 1);
        int visible = config.getInt("visible", 1);

        // Get charge configuration and check for conflicts immediately
        bool charged_only = config.getBool("charged_only", false);
        bool neutral_only = config.getBool("neutral_only", false);

        if (charged_only && neutral_only) {
            throw std::runtime_error("Configurations 'charged_only' and 'neutral_only' both set to 'YES'");
        }
        
        // Configuration option for energy distribution
        bool use_gaussian_energy = config.getBool("use_gaussian_energy", false);
        double energy_stddev = config.getDouble("energy_stddev", 0.1); // Default standard deviation
        
        // Check if energy per particle is sufficient
        if (energy <= m) {
            throw std::runtime_error("Energy (" + std::to_string(energy) + 
                                   " GeV) must be greater than particle mass (" + 
                                   std::to_string(m) + " GeV).");
        }
        
        // Get Thomson-specific parameters with updated defaults
        int thomson_iterations = config.getInt("thomson_iterations", 500);
        double thomson_step_size = config.getDouble("thomson_step_size", 0.1);

        if (N <= 0) {
            throw std::runtime_error("Number of vectors must be positive");
        }

        // Create output directory
        if (!createDirectory("./MockPileup")) {
            throw std::runtime_error("Failed to create directory: ./MockPileup");
        }

        // Construct output filename with appropriate prefix
        std::string prefix;
        if (use_gaussian_energy) {
            // Format with one decimal place for mu and sigma in the filename
            char grf_prefix[100];
            snprintf(grf_prefix, sizeof(grf_prefix), "GRF_mu%.1f_sigma%.1f_", energy, energy_stddev);
            prefix = grf_prefix;
        } else {
            prefix = "ISO_";
        }
        std::string filename = "./MockPileup/" + prefix + method + "_" + std::to_string(N) + ".dat";

        // Random number generator for charge assignment and energy
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> charge_dist(0.0, 1.0);
        
        // Gaussian distribution for energy (if enabled)
        std::normal_distribution<double> energy_dist(energy, energy_stddev);

        // Create vector generator
        auto vectorGen = createGenerator(method, gen, thomson_iterations, thomson_step_size);

        // Open output file
        std::ofstream outFile(filename);
        if (!outFile) {
            throw std::runtime_error("Failed to open output file: " + filename);
        }

        // Print headers
        outFile << std::setw(10) << "Particle"
                << std::setw(11) << "ID"
                << std::setw(11) << "Status"
                << std::setw(26) << "E"
                << std::setw(26) << "px"
                << std::setw(26) << "py"
                << std::setw(26) << "pz"
                << std::setw(26) << "m"
                << std::setw(16) << "Charge"
                << std::setw(11) << "Visible"
                << std::endl;

        // Generate vectors
        for (int i = 0; i < N; i++) {
            double dx, dy, dz;
            vectorGen->generateDirection(i, N, dx, dy, dz);
            
            // Determine energy for this particle
            double particle_energy;
            if (use_gaussian_energy) {
                // Keep drawing until we get valid energy
                do {
                    particle_energy = energy_dist(gen);
                } while (particle_energy <= m); // Ensure energy > mass
            } else {
                particle_energy = energy; // Fixed energy
            }
            
            // Calculate momentum magnitude for this energy
            double particle_p_mag = calculateMomentumMagnitude(particle_energy, m);
            
            // Calculate momentum components using momentum magnitude
            double px = particle_p_mag * dx;
            double py = particle_p_mag * dy;
            double pz = particle_p_mag * dz;
            
            // Verify mass-shell condition
            double calculatedE = std::sqrt(px*px + py*py + pz*pz + m*m);
            double energyDiff = std::abs(calculatedE - particle_energy);
            if (energyDiff > 1e-10) {  // Allow for small numerical errors
                throw std::runtime_error("Mass-shell condition E^2 - p^2 = m^2 not satisfied");
            }

            // Determine charge based on configuration
            int charge;
            if (charged_only) {
                charge = base_charge;
            } else if (neutral_only) {
                charge = 0;
            } else {
                // Randomly assign charge if neither flag is set
                charge = (charge_dist(gen) < 0.5) ? 0 : base_charge;
            }

            char buffer[512];
            snprintf(buffer, sizeof(buffer), 
                    "%10d %10d %10d %25.16e %25.16e %25.16e %25.16e %25.16e %15d %10d\n",
                    i+1, ID, status, particle_energy, px, py, pz, m, charge, visible);
            outFile << buffer;
        }

        outFile.close();
        
        // Print generation information
        std::string energy_desc = use_gaussian_energy
                              ? "Gaussian with mean " + std::to_string(energy) + " GeV and stddev " + std::to_string(energy_stddev) + " GeV"
                              : std::to_string(energy) + " GeV (fixed)";
                              
        std::cout << "Generated " << N << " vectors using " << vectorGen->getName();
        if (method == "thomson") {
            std::cout << " (iterations=" << thomson_iterations 
                     << ", step_size=" << thomson_step_size << ")";
        }
        std::cout << " method in file: " << filename << std::endl;
        std::cout << "Energy distribution: " << energy_desc << std::endl;
        if (!use_gaussian_energy) {
            std::cout << "Momentum magnitude per particle: " << calculateMomentumMagnitude(energy, m) << " GeV" << std::endl;
        }
        std::cout << "Particle mass: " << m << " GeV" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
