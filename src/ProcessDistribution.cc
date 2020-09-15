#include "ProcessDistribution.hh"

CubicSpline prob_spline;

void readDistribution() {
    std::ifstream in("distribution_E5.dat");

    std::vector<double> energy_distribution_vec_;
    std::vector<double> pdf_distribution_vec_;
    std::vector<double> cdf_distribution_vec_;

    double energy_value, pdf_value, cdf_value;
    while(in >> energy_value >> pdf_value >> cdf_value) {
        energy_distribution_vec_.push_back(energy_value);
        pdf_distribution_vec_.push_back(pdf_value);
        cdf_distribution_vec_.push_back(cdf_value);
    }

    prob_spline.SetPoints(cdf_distribution_vec_, energy_distribution_vec_);
}

double getEnergy(double randNum) {
    return prob_spline(randNum);
}