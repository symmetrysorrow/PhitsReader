#pragma once

#include <vector>
#include <exception>
#include <algorithm>
#include <Eigen/Dense>
#include<besselfi.h>
#include<iostream>
#include <fstream>
#include <sstream>
#include <numeric>  
#include "coder_array.h"
#include "struct.h"

typedef std::vector<int> vectori;
typedef std::vector<double> vectord;

void add_index_range(vectori& indices, int beg, int end, int inc = 1);

void add_index_const(vectori& indices, int value, size_t numel);

void append_vector(vectord& vec, const vectord& tail);

vectord subvector_reverse(const vectord& vec, int idx_end, int idx_start);

inline int max_val(const vectori& vec);

void filter(vectord B, vectord A, const vectord& X, vectord& Y, vectord& Zi);

void filtfilt(vectord B, vectord A, const vectord& X, vectord& Y);

std::tuple<double, double, int> peak_c(const std::vector<double>& data, int presamples, int w_max, int x_av, int w_av);

std::tuple<double, int, int> risetime(const std::vector<double>& data, double peak, int peak_index, double rise_high, double rise_low, double rate);

// Eigen::VectorXd Ç coder::array<double, 2U> Ç…ïœä∑Ç∑ÇÈä÷êî
void eigen_to_coder_array(const Eigen::VectorXd& eigen_vec, coder::array<double, 2U>& coder_arr);

std::pair<std::vector<double>, std::vector<double>> MakeCoeff(const InputParameters& InputPara);

std::tuple<int, double, double> MakeCSV(const Eigen::VectorXd& Vector, const int& Event, const std::vector<double>& b_coeff, const std::vector<double>& a_coeff);

