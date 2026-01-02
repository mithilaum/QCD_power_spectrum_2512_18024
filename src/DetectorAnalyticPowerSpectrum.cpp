#include <QtCore/QSettings>
#include <QtCore/QString>
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <limits>
#include <stdexcept>
#include <functional>

// Convert pseudorapidity to polar angle (in radians)
double etaToTheta(double eta) {
    return 2.0 * atan(exp(-eta));
}

// Improved Legendre polynomial calculation using scaled arithmetic
double legendre_improved(int l, double x) {
    if (std::abs(x) > 1.0) {
        throw std::domain_error("Argument x must be in [-1,1]");
    }
    
    if (l == 0) return 1.0;
    if (l == 1) return x;
    
    // Use scaled recursion to avoid overflow/underflow
    double p_prev2 = 1.0;
    double p_prev1 = x;
    double p_l = 0.0;
    
    const double scale = 1e-128;     // Scaling factor to prevent underflow
    const double inv_scale = 1.0/scale;
    
    for (int n = 2; n <= l; n++) {
        // Apply scaling if values getting too large or small
        if (std::abs(p_prev1) > 1e128) {
            p_prev1 *= scale;
            p_prev2 *= scale;
        } else if (std::abs(p_prev1) < 1e-128) {
            p_prev1 *= inv_scale;
            p_prev2 *= inv_scale;
        }
        
        p_l = ((2.0 * n - 1.0) * x * p_prev1 - (n - 1.0) * p_prev2) / n;
        p_prev2 = p_prev1;
        p_prev1 = p_l;
    }
    
    return p_l;
}

// Improved numerical integration using adaptive Simpson's rule
class IntegrandFunction {
private:
    int l;
public:
    IntegrandFunction(int l_value) : l(l_value) {}
    double operator()(double x) const {
        return legendre_improved(l, x);
    }
};

// Helper function for Simpson's rule integration over an interval
double simpson_rule(const IntegrandFunction& f, double a, double b) {
    double h = (b - a) / 2.0;
    double fa = f(a);
    double fm = f(a + h);
    double fb = f(b);
    return (h / 3.0) * (fa + 4.0 * fm + fb);
}

// Recursive adaptive Simpson's rule implementation with depth check
double adaptive_simpson_recursive(const IntegrandFunction& f, 
                                double a, double b, 
                                double tol, 
                                double whole,
                                int depth = 0) {
    const int max_depth = 50;  // Prevent excessive recursion
    if (depth > max_depth) {
        throw std::runtime_error("Maximum recursion depth exceeded in adaptive Simpson integration");
    }
    
    double h = (b - a) / 2.0;
    double m = a + h;
    
    double left = simpson_rule(f, a, m);
    double right = simpson_rule(f, m, b);
    double whole_new = left + right;
    
    if (std::abs(whole_new - whole) <= 15.0 * tol) {
        return whole_new + (whole_new - whole) / 15.0;
    }
    
    return adaptive_simpson_recursive(f, a, m, tol/2.0, left, depth + 1) +
           adaptive_simpson_recursive(f, m, b, tol/2.0, right, depth + 1);
}

double adaptive_simpson(int l, double a, double b, double tol = 1e-12) {
    IntegrandFunction f(l);
    double whole = simpson_rule(f, a, b);
    return adaptive_simpson_recursive(f, a, b, tol, whole, 0);
}

// Updated H_l calculation with new formulas
double calculate_Hl_improved(int l, double lambda) {
    if (lambda < 0.0 || lambda > M_PI) {
        throw std::domain_error("Lambda must be in [0,π]");
    }
    
    // For H_0: cos^2(lambda)
    if (l == 0) {
        return cos(lambda) * cos(lambda);
    }
    
    // For H_odd: 0
    if (l % 2 == 1) {
        return 0.0;
    }
    
    // For H_even (l ≥ 2): (integrate from cos(lambda) to 1 of P_l(x) dx)^2
    try {
        double integral = adaptive_simpson(l, cos(lambda), 1.0);
        return integral * integral;
    } catch (const std::exception& e) {
        std::cerr << "Error in H_l calculation for l=" << l << ": " << e.what() << std::endl;
        throw;
    }
}

// Find next non-zero H_l value with improved error checking
double find_next_nonzero(const std::vector<double>& H_values, int current_idx) {
    if (current_idx < 0 || current_idx >= static_cast<int>(H_values.size())) {
        throw std::out_of_range("Invalid current index");
    }
    
    for (int idx = current_idx + 1; idx < static_cast<int>(H_values.size()); ++idx) {
        if (H_values[idx] > std::numeric_limits<double>::epsilon()) {
            return H_values[idx];
        }
    }
    return 0.0;
}

int main(int argc, char* argv[]) {
    try {
        if (argc != 2) {
            throw std::runtime_error("Usage: " + std::string(argv[0]) + " <config_file>");
        }

        // Parse configuration file using QSettings
        QSettings config(QString::fromStdString(argv[1]), QSettings::IniFormat);
        if (config.status() != QSettings::NoError) {
            throw std::runtime_error("Error reading configuration file");
        }
        
        // Read parameters with validation
        int l_max = config.value("power/lMax", 256).toInt();
        if (l_max < 0) {
            throw std::runtime_error("l_max must be non-negative");
        }
        
        double eta_max = config.value("detector/etaMax_cal", 2.0).toDouble();
        if (eta_max <= 0.0) {
            throw std::runtime_error("eta_max must be positive");
        }
        
        // Read normalization preference
        QString normalizeStr = config.value("power/normalizeToH0", "NO").toString().toUpper();
        bool normalize = (normalizeStr == "YES");
        
        // Convert eta to theta
        double theta = etaToTheta(eta_max);
        double lambda = theta;

        // Create output filename
        std::stringstream ss;
        ss << "powerspectrum_detector_holes_etaMax_" << std::fixed << std::setprecision(3) 
           << eta_max << (normalize ? "_normalized" : "_raw") << ".dat";
        std::string outfile_name = ss.str();

        // Open output file
        std::ofstream outfile(outfile_name);
        if (!outfile.is_open()) {
            throw std::runtime_error("Could not create output file: " + outfile_name);
        }

        // Set output precision to maximum available for double
        outfile.precision(std::numeric_limits<double>::max_digits10);
        outfile.setf(std::ios::scientific);

        // Calculate H_0 first if normalization is needed
        double H_0 = 1.0;
        if (normalize) {
            H_0 = calculate_Hl_improved(0, lambda);
            if (H_0 <= std::numeric_limits<double>::epsilon()) {
                throw std::runtime_error("H_0 is zero or negative, cannot normalize");
            }
        }

        // Calculate all H_l values
        std::vector<double> H_values(l_max + 1);
        const int progress_interval = std::max(1, l_max / 20);  // Show progress every 5%
        
        for (int l = 0; l <= l_max; ++l) {
            try {
                H_values[l] = calculate_Hl_improved(l, lambda);
                if (normalize) {
                    H_values[l] /= H_0;
                }
                
                // Show progress
                if (l % progress_interval == 0) {
                    std::cout << "Computing l = " << l << " / " << l_max << "\n";
                }
            } catch (const std::exception& e) {
                throw std::runtime_error("Error calculating H_" + std::to_string(l) + ": " + e.what());
            }
        }

        // Output values with running averages
        for (int l = 0; l <= l_max; ++l) {
            double current = H_values[l];
            double next = find_next_nonzero(H_values, l);
            double running_avg = (current + next) / 2.0;
            outfile << l << " " << current << " " << running_avg << "\n";
        }
        
        outfile.close();
        std::cout << "\nResults written to " << outfile_name << "\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
