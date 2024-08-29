#pragma once
#include <cmath>
#include<complex>
#include <iostream>
#include<Eigen/Dense>
#include <iomanip>

//intと粒子の種類の対応
std::map<int, std::string> itype = {
    {12, "electron"},
    {13, "positron"},
    {14, "photon"}
};
//数字を受け取り対応する粒子を返す関数
inline std::string GetItype(const double& ityp) {
    auto it = itype.find(static_cast<int>(ityp));
    if (it != itype.end())
    {
        return it->second;
    }
    else
    {
        return "unknown";  // デフォルト値を返す
    }
};
//Eventに関する構造体
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

//空白で文章を分割する関数
std::vector<double> split_line(const std::string& line) {
    std::vector<double> column;
    std::istringstream stream(line);
    std::string token;

    while (stream >> token) {  // 空白をスキップしてトークンを取得
        try {
            column.push_back(std::stof(token));  // トークンをdoubleに変換
        }
        catch (const std::invalid_argument& e) {
            std::cerr << "Invalid argument: " << token << " cannot be converted to double." << std::endl;
        }
        catch (const std::out_of_range& e) {
            std::cerr << "Out of range: " << token << " is out of range for double." << std::endl;
        }
    }

    return column;
}

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

InputParameters ReadInputJson(const std::string& InputPath){

    InputParameters InputPara;
    //jsonの読み込み
    std::ifstream InputStream(InputPath);

    if (!InputStream.is_open()) {
        std::cerr << "ファイルを開くことができません: " << InputPath << std::endl;
    }

    nlohmann::json InputJson;
    try {
        InputJson = nlohmann::json::parse(InputStream);
    }
    catch (const nlohmann::json::parse_error& e) {
        std::cerr << "JSONの解析中にエラーが発生しました: " << e.what() << std::endl;
    }

    InputPara.C_abs = InputJson["C_abs"];
    InputPara.C_tes = InputJson["C_tes"];
    InputPara.G_abs_abs = InputJson["G_abs-abs"];
    InputPara.G_abs_tes = InputJson["G_abs-tes"];
    InputPara.G_tes_bath = InputJson["G_tes-bath"];
    InputPara.R = InputJson["R"];
    InputPara.R_l = InputJson["R_l"];
    InputPara.T_c = InputJson["T_c"];
    InputPara.T_bath = InputJson["T_bath"];
    InputPara.alpha = InputJson["alpha"];
    InputPara.beta = InputJson["beta"];
    InputPara.L = InputJson["L"];
    InputPara.n = InputJson["n"];
    InputPara.E = InputJson["E"];
    InputPara.length = InputJson["length"];
    InputPara.n_abs = InputJson["n_abs"];
    InputPara.rate = InputJson["rate"];
    InputPara.samples = InputJson["samples"];
    InputPara.data_samples = InputJson["data_samples"];
    InputPara.cutoff = InputJson["cutoff"];
    InputPara.history = InputJson["history"];
    InputPara.output = InputJson["output"];

    return InputPara;
}

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

int InBlock(const std::vector<double>& Block, const double& x_deposit, const double& y_deposit, const double& z_deposit)
{
    // std::cout << "x_deposit:" << x_deposit << std::setprecision(8) << "\ny_deposit:" << y_deposit << std::setprecision(8) << "\nz_deposit:" << z_deposit << std::setprecision(8) << "\n";
     //std::cout << "block:[";
     //std::cout << "]\n";
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

