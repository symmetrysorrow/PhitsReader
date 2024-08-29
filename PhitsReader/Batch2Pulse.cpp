#pragma once

#include "Batch2Pulse.h"
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
