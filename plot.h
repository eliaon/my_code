#pragma once
#ifndef PLOT_H
#define PLOT_H

#include <string>

std::string get_meson();

void read_sigma_exp(
    const std::string& filename,
    std::vector<int>& dataset,
    std::vector<double>& W,
    std::vector<double>& sigma,
    std::vector<double>& err);

void plot_dsigma_dt(std::string meson);

void plot_N_models(const std::vector<std::pair<std::string,std::string>>& files, double x);

void compare_N_models(double x);

void plot_overlap();

void plot_rapidity();

void plot_sigma(std::string meson, std::string csv_file, std::string dipolemodel = "GBW");

void plot_sigma_phi(std::string csv_file, std::string dipolemodel);
void plot_sigma_Jpsi(std::string csv_file, std::string dipolemodel);

void plot_sigma_Jpsi_gammaPb_GBW(std::string csv_file);
void plot_sigma_phi_gammaPb_GBW(std::string csv_file);

void plot_TA_b();

void plot_rapidez_PbPb_Jpsi(std::string csv_file, double sqrt_s);
void plot_rapidez_PbPb_phi(std::string csv_file, double sqrt_s);

void plot_N_dglap(std::string csv_file1, std::string csv_file2);

void plot_overlap_fc(std::string csv_file, std::string csv_file_fc, std::string meson);

void plot_xf_multiQ2();
void plot_N_multi(
    const std::vector<std::string>& filenames,
    const std::vector<double>& x_vals,
    const std::string& modelo
);

#endif // PLOT_H
