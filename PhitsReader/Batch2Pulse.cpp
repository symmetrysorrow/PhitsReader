#pragma once

#include "Batch2Pulse.h"
#include<filesystem>
#include <iostream>

#include"Dump2Batch.h"

PulseParameters::PulseParameters(const InputParameters& InputPara)
    : C_abs(InputPara.C_abs / InputPara.n_abs),
    G_abs_abs(InputPara.G_abs_abs* (InputPara.n_abs - 1)),
    I(std::sqrt((InputPara.G_tes_bath* InputPara.T_c* (1 - std::pow(InputPara.T_bath / InputPara.T_c, InputPara.n))) / (InputPara.n * InputPara.R))),
    t_el(InputPara.L / (InputPara.R_l + InputPara.R * (1 + InputPara.beta))),
    L_I((InputPara.alpha* std::pow(I, 2)* InputPara.R) / (InputPara.G_tes_bath * InputPara.T_c)),
    t_I(InputPara.C_tes / ((1 - L_I) * InputPara.G_tes_bath)) {}

// 等差数列を作成する関数
std::vector<double> linspace(double start, double stop, int num) {
    std::vector<double> values;
    if (num <= 0) return values;


    double step = (stop - start) / (num - 1);

    for (int i = 0; i < num; ++i) {
        values.push_back(start + i * step);
    }

    return values;
}

// なぞ関数1
int InBlock(const std::vector<double>& Block, const double& x_deposit, const double& y_deposit, const double& z_deposit) {
    if ((-1.0 <= x_deposit || x_deposit <= 1) && (-0.05 <= y_deposit || y_deposit <= 0.05) && (-0.1 <= z_deposit || z_deposit <= 0.0)) {
        for (int i = 0; i < Block.size(); ++i) {
            if (Block[i] >= x_deposit) {
                return i;
            }
        }
    }
    return 0;
}

Eigen::MatrixXd MakeMatrix_M(const PulseParameters& PulsePara, const InputParameters& InputPara)
{
    const int n_abs = InputPara.n_abs;
	const int n_abs_1 = n_abs + 1;
	const int n_abs_2 = n_abs + 2;
	const int n_abs_3 = n_abs + 3;
	const int n_abs_4 = n_abs + 4;


	Eigen::MatrixXd Matrix_M(n_abs_4, n_abs_4);
	Matrix_M = Eigen::MatrixXd::Zero(n_abs_4, n_abs_4);

	for (int i = 0; i < n_abs_4; i++)
	{
		for (int j = 0; j < n_abs_4; j++)
		{
			if (j == i - 1)
			{
				Matrix_M(i, j) = -PulsePara.G_abs_abs / PulsePara.C_abs;
			}
			if (j == i)
			{
				Matrix_M(i, j) = 2 * PulsePara.G_abs_abs / PulsePara.C_abs;
			}
			if (j == i + 1)
			{
				Matrix_M(i, j) = -PulsePara.G_abs_abs / PulsePara.C_abs;
			}
		}
	}

	Matrix_M(0, 0) = 1 / PulsePara.t_el;
	Matrix_M(0, 1) = PulsePara.L_I * InputPara.G_tes_bath / (PulsePara.I * InputPara.L);

	Matrix_M(1, 0) = -PulsePara.I * InputPara.R * (2 + InputPara.beta) / InputPara.C_tes;
	Matrix_M(1, 1) = 1 / PulsePara.t_I + (InputPara.G_abs_tes / InputPara.C_tes);
	Matrix_M(1, 2) = -InputPara.G_abs_tes / InputPara.C_tes;

	Matrix_M(2, 1) = -InputPara.G_abs_tes / PulsePara.C_abs;
	Matrix_M(2, 2) = InputPara.G_abs_tes / PulsePara.C_abs + PulsePara.G_abs_abs / PulsePara.C_abs;
	Matrix_M(2, 3) = -PulsePara.G_abs_abs / PulsePara.C_abs;

	Matrix_M(n_abs_1, n_abs) = -PulsePara.G_abs_abs / PulsePara.C_abs;
	Matrix_M(n_abs_1, n_abs_1) = InputPara.G_abs_tes / PulsePara.C_abs + PulsePara.G_abs_abs / PulsePara.C_abs;
	Matrix_M(n_abs_1, n_abs_2) = -InputPara.G_abs_tes / PulsePara.C_abs;

	Matrix_M(n_abs_2, n_abs_1) = -InputPara.G_abs_tes / InputPara.C_tes;
	Matrix_M(n_abs_2, n_abs_2) = 1 / PulsePara.t_I + InputPara.G_abs_tes / InputPara.C_tes;
	Matrix_M(n_abs_2, n_abs_3) = -PulsePara.I * InputPara.R * (2 + InputPara.beta) / InputPara.C_tes;

	Matrix_M(n_abs_3, n_abs_2) = PulsePara.L_I * InputPara.G_tes_bath / (PulsePara.I * InputPara.L);
	Matrix_M(n_abs_3, n_abs_3) = 1 / PulsePara.t_el;

	Matrix_M *= -1;

	return Matrix_M;
}

Eigen::MatrixXd MakeMatrix_X(const PulseParameters& PulsePara, const InputParameters& InputPara, const std::vector<int>& pixel)
{
	const int n_abs = InputPara.n_abs;
	const int n_abs_1 = n_abs + 1;
	const int n_abs_2 = n_abs + 2;
	const int n_abs_3 = n_abs + 3;
	const int n_abs_4 = n_abs + 4;

	Eigen::MatrixXd Matrix_X(n_abs_2, n_abs_4);
	Matrix_X = Eigen::MatrixXd::Zero(n_abs_2, n_abs_4);

	const double Same = InputPara.E * 1e3 * PulsePara.e_const / InputPara.C_tes;

	for (const int& pix : pixel)
	{
		if (pix <= n_abs_2)
		{
			Matrix_X(pix, pix + 1) = InputPara.positions[pix - 1] * PulsePara.e_const / PulsePara.C_abs;
		}
	}

	Matrix_X(0, 1) = Same;
	Matrix_X(n_abs_1, n_abs_2) = Same;

	return Matrix_X;
}
