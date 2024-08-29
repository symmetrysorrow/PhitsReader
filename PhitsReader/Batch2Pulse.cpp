#pragma once
#include <vector>
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