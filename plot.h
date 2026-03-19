#pragma once
#ifndef PLOT_H
#define PLOT_H

#include <string>

std::string get_meson();

void plot_dsigma_dt(std::string meson);

void plot_N_models();

void plot_overlap();

void plot_rapidity();

void plot_sigma(std::string meson, std::string csv_file);

void plot_sigma_phi(std::string csv_file);
void plot_sigma_Jpsi(std::string csv_file);

#endif // PLOT_H
