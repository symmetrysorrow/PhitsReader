#pragma once
#include <Eigen/Dense>

Eigen::VectorXd  inverseFFTNoise(const Eigen::VectorXcd& noise_spe, int noise_samples, int samples);

Eigen::VectorXcd random_noise(const Eigen::VectorXd& spe, unsigned seed);

void AddNoise(const Eigen::VectorXd& Noise_dense, Eigen::VectorXd& Pulse);

Eigen::VectorXd readLinesToEigen(const std::string& filePath);