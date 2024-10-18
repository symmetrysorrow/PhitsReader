#pragma once
#include<vector>

struct EventInfo {
    int ityp;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<double> E;
    std::vector<double> x_deposit;
    std::vector<double> y_deposit;
    std::vector<double> z_deposit;
    std::vector<double> E_deposit;
};

// Input.jsonの各パラメータの構造体
struct InputParameters
#if 1
{
    double C_abs;
    double C_tes;
    double G_abs_abs;
    double G_abs_tes;
    double G_tes_bath;
    double R;
    double R_l;
    double T_c;
    double T_bath;
    double alpha;
    double beta;
    double L;
    double n;
    double E;
    double length;
    double height;
    double depth;
    int n_abs;
    double rate;
    double samples;
    std::vector<int> positions;
    int data_samples;
    int cutoff;
    int history;
    double output;
    bool noise;
    bool SavePulse;
};
#endif