// Batch2Pulse.hpp などのヘッダーファイル
#ifndef BATCH2PULSE_HPP
#define BATCH2PULSE_HPP

#pragma once
#include <vector>
//#include"Dump2Batch.hpp"
#include<Eigen/Dense>

struct InputParameters;

//パルスに関する構造体
struct PulseParameters {
    const double k_b = 1.381 * 1.0e-23;
    const double ptfn_Flink = 0.5;
    const double e_const = 1.602e-19;

    const double C_abs;
    const double G_abs_abs;
    const double I;
    const double t_el;
    const double L_I;
    const double t_I;

    PulseParameters(const InputParameters& InputPara):
	C_abs (InputPara.C_abs / InputPara.n_abs),
	G_abs_abs(InputPara.G_abs_abs * (InputPara.n_abs - 1)),
    I(std::sqrt((InputPara.G_tes_bath * InputPara.T_c * (1 - std::pow(InputPara.T_bath / InputPara.T_c, InputPara.n))) / (InputPara.n * InputPara.R))),
    t_el(InputPara.L / (InputPara.R_l + InputPara.R * (1 + InputPara.beta))),
    L_I((InputPara.alpha * std::pow(I, 2) * InputPara.R) / (InputPara.G_tes_bath * InputPara.T_c)),
    t_I(InputPara.C_tes / ((1 - L_I) * InputPara.G_tes_bath)){}
    
};

//等差数列を作成する関数
std::vector<double> linspace(double start, double stop, int num) {
    std::vector<double> values;
    if (num <= 0) return values; // numが0以下の場合は空の配列を返す

    // ステップ幅を計算
    double step = (stop - start) / (num - 1);

    // 値を生成
    for (int i = 0; i < num; ++i) {
        values.push_back(start + i * step);
    }

    return values;
}
//なぞ関数1
int InBlock(const std::vector<double>& Block, const double& x_deposit, const double& y_deposit, const double& z_deposit)
{
    if ((-1.0 <= x_deposit || x_deposit <= 1) && (-0.05 <= y_deposit || y_deposit <= 0.05) && (-0.1 <= z_deposit || z_deposit <= 0.0))
    {
        for (int i = 0; i < Block.size(); ++i) {
            if (Block[i] >= x_deposit) {
                return i; // 条件を満たす最初の要素のインデックス
            }
        }
    }
    return 0;
}

#endif