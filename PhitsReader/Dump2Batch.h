#pragma once

#include <filesystem>
#include <map>
#include <vector>

// intと粒子の種類の対応マップ
extern std::map<int, std::string> itype;

// 数字を受け取り対応する粒子を返す関数
std::string GetItype(const double& ityp);

// Eventに関する構造体
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

// 空白で文章を分割する関数
std::vector<double> split_line(const std::string& line);

// Input.jsonの各パラメータの構造体
struct InputParameters
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
    int n_abs;
    double rate;
    double samples;
    std::vector<double> positions;
    int data_samples;
    int cutoff;
    int history;
    std::string output;
};

// Input.jsonを読み出す関数
InputParameters ReadInputJson(const std::string& InputPath);

// dumpall.datをbatchにする関数
std::map<int, std::map<int, EventInfo>> ReadDump(const std::string& DumpPath);

// batchをoutput.jsonに書き出す関数
void WriteOutput(const std::map<int, std::map<int, EventInfo>>& batch, const std::string& output_file);
