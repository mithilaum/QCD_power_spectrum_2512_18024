#include "Post_PYTHIA_PowerJets.hpp"
#include "Admin.hpp"
#include "kdpVectors.hpp"
#include "PowerSpectrum.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <QtCore/QSettings>

// ------------------------------------------------------------
// ParticleData + helpers
// ------------------------------------------------------------
struct ParticleData {
    kdp::Vector4<double> momentum;
    bool isCharged;
    bool isVisible;
};

std::vector<ParticleData> readParticlesFromFile(const std::string& filename) {
    std::vector<ParticleData> particles;
    std::ifstream file(filename);
    if (!file.is_open()) return particles;

    std::string line;
    std::getline(file, line); // header

    while (std::getline(file, line)) {
        if (line.empty() || line.find('#') != std::string::npos) continue;

        std::istringstream iss(line);
        int particle_num, id, status;
        double E, px, py, pz, mass, charge, visible;

        if (iss >> particle_num >> id >> status
                >> E >> px >> py >> pz >> mass
                >> charge >> visible) {

            ParticleData p;
            kdp::Vector3<double> p3(px, py, pz);
            p.momentum  = kdp::Vector4<double>(E, p3, kdp::Vec4from2::Energy);
            p.isCharged = (charge != 0.0);
            p.isVisible = (visible > 0.5);
            particles.push_back(p);
        }
    }
    return particles;
}

// ------------------------------------------------------------
// Helper: filter visible particles + beam holes
// ------------------------------------------------------------
std::vector<kdp::Vector4<double>>
filterVisibleParticles(const std::vector<ParticleData>& data,
                       bool useBeamHoles,
                       double etaMax)
{
    std::vector<kdp::Vector4<double>> out;
    for (auto const& p : data) {
        if (!p.isVisible) continue;
        if (useBeamHoles &&
            std::fabs(p.momentum.p().Eta()) > etaMax)
            continue;
        out.push_back(p.momentum);
    }
    return out;
}

// ------------------------------------------------------------
// Compute <f|f> for particles entering power spectrum
// ------------------------------------------------------------
double compute_ff(const std::vector<kdp::Vector4<double>>& particles) {
    double ff = 0.0;
    auto phat = PowerSpectrum::PhatF::To_PhatF_Vec(particles);
    for (auto const& p : phat) ff += p.f * p.f;
    return ff;
}

// ------------------------------------------------------------
// main()
// ------------------------------------------------------------
int main() {

    const std::string configFilePath =
        "/home/mithila/EmPowerSpectrum/CONF/EvenOdd_PowerSpectrum.conf";

    QSettings cfg(configFilePath.c_str(), QSettings::NativeFormat);

    std::string file_ME =
        cfg.value("main/inputDataFile_hard_process_final", "")
           .toString().toStdString();

    std::string file_all =
        cfg.value("main/inputDataFile_all", "")
           .toString().toStdString();

    if (file_ME.empty() || file_all.empty()) {
        std::cerr << "Both inputDataFile_hard_process_final and "
                     "inputDataFile_all must be specified\n";
        return EXIT_FAILURE;
    }

    bool useBeamHoles =
        cfg.value("detector/beam_holes", "ON")
           .toString().toUpper() == "ON";

    double etaMax =
        cfg.value("detector/etaMax_cal", 4.9).toDouble();

    // ------------------------------------------------------------
    // Read particles
    // ------------------------------------------------------------
    auto particles_ME  = readParticlesFromFile(file_ME);
    auto particles_all = readParticlesFromFile(file_all);

    if (particles_ME.empty() || particles_all.empty()) {
        std::cerr << "Failed to read one of the input files\n";
        return EXIT_FAILURE;
    }

    // ------------------------------------------------------------
    // Initialize PowerJets
    // ------------------------------------------------------------
    Post_PYTHIA_PowerJets pj_ME(configFilePath);
    Post_PYTHIA_PowerJets pj_all(configFilePath);

    if (pj_ME.GetStatus() != Post_PYTHIA_PowerJets::Status::OK ||
        pj_all.GetStatus() != Post_PYTHIA_PowerJets::Status::OK) {
        std::cerr << "PowerJets initialization failed\n";
        return EXIT_FAILURE;
    }

    auto filtered_ME  = filterVisibleParticles(particles_ME,  useBeamHoles, etaMax);
    auto filtered_all = filterVisibleParticles(particles_all, useBeamHoles, etaMax);

    pj_ME.SetParticles(filtered_ME);
    pj_all.SetParticles(filtered_all);

    // ------------------------------------------------------------
    // Retrieve spectra
    // ------------------------------------------------------------
    const auto& Hl_ME = pj_ME.Get_Hl_FinalState();
    const auto& Hl_FS = pj_all.Get_Hl_FinalState();

    size_t maxL = std::min(Hl_ME.size(), Hl_FS.size());
    maxL -= maxL % 2;

    // ------------------------------------------------------------
    // Output paths
    // ------------------------------------------------------------
    std::string outputDir =
        cfg.value("main/outputDir", "./").toString().toStdString();
    if (outputDir.back() != '/') outputDir += "/";

    std::string tag =
        useBeamHoles ? "_with_holes_" + std::to_string(etaMax)
                     : "_spherical";

    std::string process = GetProcess(file_all);

    std::string outData =
        outputDir + "H_" + process + tag + "_evenOdd.dat";

    std::string outLog =
        outputDir + "log_" + process + tag + "_evenOdd.txt";

    // ------------------------------------------------------------
    // Write EVEN/ODD DATA FILE
    // ------------------------------------------------------------
    std::ofstream out(outData, std::ios::trunc);
    out << "# l  H_ME_l   H_ME_(l+1)   H_finalState_l   H_finalState_(l+1)\n";
    out << "  -1 nan 1. nan 1.\n";

    out << std::scientific << std::setprecision(16);

    for (size_t l = 0; l < maxL; l += 2) {
        out << std::setw(4) << (l + 1) << "  "
            << Hl_ME[l]  << "  "
            << Hl_ME[l+1] << "  "
            << Hl_FS[l]  << "  "
            << Hl_FS[l+1] << "\n";
    }
    out.close();

    // ------------------------------------------------------------
    // Write LOG FILE
    // ------------------------------------------------------------
    double ff_ME  = compute_ff(filtered_ME);
    double ff_FS  = compute_ff(filtered_all);

    std::ofstream log(outLog, std::ios::trunc);
    if (log.is_open()) {

      log << std::scientific << std::setprecision(16);
     
      log << "Summary:\n";
      log << "Input file (ME): " << file_ME << "\n";
      log << "Input file (final state): " << file_all << "\n";
      log << "Output file: " << outData << "\n\n";
      
      log << "Number of particles (ME): "
	  << filtered_ME.size() << "\n";
      log << "Partons <f|f>: "
	  << ff_ME << "\n\n";
      
      log << "Number of particles (final state): "
	  << filtered_all.size() << "\n";
      log << "Particles <f|f>: "
	  << ff_FS << "\n\n";
      
      log << "Maximum l value: " << maxL << "\n";
      log.close();
    }
    
    std::cout << "Wrote even/odd spectrum to:\n  " << outData << "\n";
    std::cout << "Log written to:\n  " << outLog << "\n";

    return 0;
}
