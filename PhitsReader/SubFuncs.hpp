#pragma once
#include <Eigen/Dense>

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

	if((1.0<=x_deposit||x_deposit<=1) && (-0.05<=y_deposit||y_deposit<=0.05)&& (-0.1<=z_deposit||z_deposit<=0.0))
	{
        for (size_t i = 0; i < Block.size(); ++i) {
            if (Block[i] >= x_deposit) {
                return i; // 条件を満たす最初の要素のインデックス
            }
        }
	}
    return 0;
}

