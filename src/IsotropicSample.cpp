//CM: Updated by Mithila. 
// Copyright (C) 2018, Keith Pedersen (Keith.David.Pedersen@gmail.com)
// All rights reserved

#include "kdpVectors.hpp"
#include "kdpStdVectorMath.hpp"
#include "pqRand.hpp"
#include "distributions.hpp"
#include <cmath>
#include <sys/stat.h>
#include <errno.h>

#include <QtCore/QSettings>
#include <QtCore/QStringList>

#include "PowerSpectrum.hpp"
#include "ArrogantDetector.hpp"

#include "kdpHistogram.hpp"

#include "../src/IsotropicVoronoi.cpp"

// Function to create directory if it doesn't exist
bool createDirectory(const std::string& path) {
#ifdef _WIN32
    return _mkdir(path.c_str()) == 0 || errno == EEXIST;
#else
    return mkdir(path.c_str(), 0777) == 0 || errno == EEXIST;
#endif
}

int main()
{
    kdp::LimitMemory(2.);
    
    QSettings parsedINI("IsotropicSample.ini", QSettings::IniFormat);
    
    ArrogantDetector_Lepton detector(parsedINI);
    
    pqRand::engine gen(false); // don't seed
    gen.Seed_Reuse("pqRand.seed");
    
    // Create output directory
    const std::string outputDir = "output_files";
    if (!createDirectory(outputDir)) {
        std::cerr << "Error: Could not create output directory: " << outputDir << std::endl;
        return 1;
    }
    
    // Create a log of important other information for the run. Start with the ini file
    std::string logFileName = parsedINI.value("logName", "default.log").toString().toStdString();
    
    {
        std::ofstream runLog(outputDir + "/" + logFileName, std::ios::trunc | std::ios::binary);
        std::ifstream ini("IsotropicVoronoi.ini", std::ios::binary);
        
        runLog << ini.rdbuf();
    }
    
    std::ofstream runLog(outputDir + "/" + logFileName, std::ios::trunc);
    runLog << "\n\n";
    runLog << "=============== INI file included above ==================\n";  
    
    /////////////////////////////////////////////////////////////////////
    
    // The number of iso particles sampled, per generated particles, to determine Voronoi area
    size_t const sampleSize_area_perParticle = 
        size_t(parsedINI.value("sampleSize_area_perParticle", 1e4).toDouble());
    
    size_t const lMax = size_t(parsedINI.value("lMax", 1024).toInt());
    
    auto sampleSizeVec = parsedINI.value("sampleSizes").toStringList();    
    
    PowerSpectrum Hl_computer;
    /////////////////////////////////////////////////////////////////////
    kdp::BinSpecs fBins("fBins", 200, {0., 5.});
    kdp::BinSpecs fBins_log("fBins_log", 200, {1e-6, 1e3});
    
    /////////////////////////////////////////////////////////////////////
    char buff[1024];    
    sprintf(buff, "squareWidth: %.3e\n", detector.GetSettings().squareWidth);
        runLog << buff;
    sprintf(buff, "etaMax (trk/twr): %.3e  %.3e\n", detector.GetSettings().etaMax_track, 
        detector.GetSettings().etaMax_cal);
        runLog << buff;
        
    // sprintf(buff, "nTowers: %lu\n\n", detector.NumTowers()); // CM: Changed by Mithila to below.
	sprintf(buff, "nTowers: %lu\n\n", detector.Towers().size());	
        runLog << buff;
    
    for(auto const& QsampleSize : sampleSizeVec)
    {
        size_t const sampleSize = size_t(QsampleSize.toInt());
        
        /////////////////////////////////////////////////////////////////////        
        kdp::Histogram1 f_hist(outputDir + "/f_" + std::to_string(sampleSize), fBins);
        kdp::Histogram1 f_hist_log(outputDir + "/f_log_" + std::to_string(sampleSize), fBins_log);
        
        sprintf(buff, "-------------------------------\n%lu\n", sampleSize);
            runLog << buff;
        
        auto const sample = MakeIsoVoronoiVectors(gen, sampleSize, 
            sampleSize_area_perParticle * sampleSize, true);
            
        for(auto const& vec : sample)
        {
            double const f = vec.Mag();
            f_hist.Fill(f*double(sampleSize));
            f_hist_log.Fill(f*double(sampleSize));
        }
            
        // Convert to massless 4-vectors to give to detector
        std::vector<kdp::Vec4> sampleP4(sample.begin(), sample.end());
        
        // Put the particles in as both tracks and towers
        // Because they are massless, there is no contamination of the calorimeter from 
        // the E - |p>| of the track-subtracted towers.
        detector.Clear();
            // detector.PartialFill(sampleP4, sampleP4); // CM: Changed by Mithila to below.
            detector(sampleP4, sampleP4);	
        // detector.Finalize(false); // CM: Commented out by Mithila since operator() handles the finalization.
        
        std::vector<std::vector<double>> Hl_set;
        
        {   
            auto track_NN = NearestNeighborDistance(detector.Tracks());
            
            double mean = 0.; // Mean from f_i
            double mean_v2 = 0.; // Mean from f_i*f_j
            double mean_geo = 0.; // Geometric mean from w_ij
            double ff = 0.;
                        
            {
                auto tracks = PowerSpectrum::PhatF::To_PhatF_Vec(detector.Tracks());
                    assert(track_NN.size() == tracks.size());
                double weight = 0.;
                
                for(size_t i = 0; i < track_NN.size(); ++i)
                {
                    double const weight_i = tracks[i].f * tracks[track_NN[i].second].f;
                    
                    mean += track_NN[i].first * tracks[i].f;
                    mean_v2 += weight_i * track_NN[i].first;
                    mean_geo += weight_i * std::log(track_NN[i].first);
                    
                    ff += kdp::Squared(tracks[i].f);
                    weight += weight_i;
                }               
                mean_v2 /= weight;
                mean_geo = std::exp(mean_geo / weight);
            }
            std::sort(track_NN.begin(), track_NN.end());
            
            double const median = track_NN[track_NN.size() / 2].first;
            
            Hl_set.push_back(Hl_computer.Hl_Obs(lMax, PowerSpectrum::ShapedParticleContainer(detector.Tracks())));
            
            Hl_set.push_back(Hl_computer.Hl_Obs(lMax, PowerSpectrum::ShapedParticleContainer(detector.Tracks(), 
                h_Gaussian(std::sin(0.25*mean)*std::sqrt(-2./std::log1p(-0.9))))));     
        
            sprintf(buff, "\n\ntracks:    %4lu\nmin NN:    %.3e  \nmedian NN: %.3e  \nmean NN (f_i*angle):   %.3e  \nmean NN (w_ij*angle): %.3e  \nmean (geometric, w_ij) NN: %.3e  \nAngularResolution NN: %.3e  \nmax NN:    %.3e  \nexpected   %.3e  \n <f|f> :   %.3e  \n", 
                detector.Tracks().size(), 
                track_NN.front().first,
                median,
                mean,
                mean_v2,
                mean_geo,
                PowerSpectrum::AngularResolution(PowerSpectrum::PhatF::To_PhatF_Vec(detector.Tracks()), 1.),
                track_NN.back().first,
                4.*std::asin(sqrt(1./double(detector.Tracks().size()))),
                ff);
            runLog << buff;
        }
            
        //////  
        std::vector<kdp::Vec3> towers(detector.Towers().begin(), detector.Towers().end());
        
        {
            auto tower_NN = NearestNeighborDistance(towers);
            
            double mean = 0.;
            double mean_v2 = 0.;
            double mean_geo = 0.;
            double ff = 0.;
                    
            {
                auto towers_p = PowerSpectrum::PhatF::To_PhatF_Vec(towers);
                double weight = 0.;
                
                assert(tower_NN.size() == towers_p.size());
                
                for(size_t i = 0; i < tower_NN.size(); ++i)
                {
                    double const weight_i = towers_p[i].f * towers_p[tower_NN[i].second].f;
                    
                    mean += tower_NN[i].first * towers_p[i].f;
                    mean_v2 += weight_i * tower_NN[i].first;
                    mean_geo += weight_i * std::log(tower_NN[i].first);
                    
                    ff += kdp::Squared(towers_p[i].f);
                    weight += weight_i;
                }
                mean_v2 /= weight;
                mean_geo = std::exp(mean_geo / weight);
            }       
            std::sort(tower_NN.begin(), tower_NN.end());
            
            double const median = tower_NN[tower_NN.size() / 2].first;
            
            Hl_set.push_back(Hl_computer.Hl_Obs(lMax, PowerSpectrum::ShapedParticleContainer(towers)));
            
            Hl_set.push_back(Hl_computer.Hl_Obs(lMax, PowerSpectrum::ShapedParticleContainer(towers, 
                h_Gaussian(std::sin(0.25*mean)*std::sqrt(-2./std::log1p(-0.9))))));
            
            sprintf(buff, "\n\ntowers:    %4lu\nmin NN:    %.3e  \nmedian NN: %.3e  \nmean NN (f_i*angle):   %.3e  \nmean NN (w_ij*angle): %.3e  \nmean (geometric, w_ij) NN: %.3e  \nAngularResolution NN: %.3e  \nmax NN:    %.3e  \nexpected   %.3e  \n <f|f> :   %.3e  \n", 
                towers.size(), 
                tower_NN.front().first, 
                median,
                mean,
                mean_v2,
                mean_geo,
                PowerSpectrum::AngularResolution(PowerSpectrum::PowerSpectrum::PhatF::To_PhatF_Vec(towers), 1.),
                tower_NN.back().first,
                4.*std::asin(sqrt(1./double(towers.size()))),
                ff);
                
            runLog << buff;
        }
                
        std::string fileName = outputDir + "/" + QsampleSize.toStdString() + "_" + 
            std::to_string(int(std::round(kdp::ToDegrees(detector.GetSettings().squareWidth)))) + ".dat";
        
        PowerSpectrum::WriteToFile(fileName, Hl_set);
        runLog.flush();
    }
    
    return 0;
}
