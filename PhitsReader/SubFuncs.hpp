#pragma once
#include <cmath>
#include<complex>
//#EIGEN_USE_MKL_ALL
#include <iostream>
#include<Eigen/Dense>

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
//パルスを作る？関数
void Model(const InputParameters& InputPara, Eigen::MatrixXd& I_t0, Eigen::MatrixXd& I_t1)
{
    constexpr double k_b = 1.381 * 1.0e-23;
    constexpr double ptfn_Flink = 0.5;
    constexpr double e = 1.602e-19;

    const double C_abs = InputPara.C_abs / InputPara.n_abs;
    const double G_abs_abs = InputPara.G_abs_abs * (InputPara.n_abs - 1);
    const double I=std::sqrt(InputPara.G_tes_bath* InputPara.T_c*(1-std::pow((InputPara.T_bath/InputPara.T_c),InputPara.n)))/(InputPara.n*InputPara.R);
    const double t_el = InputPara.L / (InputPara.R_l + InputPara.R * (1 + InputPara.beta));
    const double L_I = (InputPara.alpha * std::pow(I, 2) * InputPara.R) / (InputPara.G_tes_bath * InputPara.T_c);
    const double t_I = InputPara.C_tes / ((1 - L_I) * InputPara.G_tes_bath);

    const int n_abs_1 = InputPara.n_abs + 1;
    const int n_abs_2 = InputPara.n_abs + 2;
    const int n_abs_3 = InputPara.n_abs + 3;
    const int n_abs_4 = InputPara.n_abs + 4;

    std::vector<int> pixel;
    int count = 1;
    for(const double& pos:InputPara.positions)
    {
	    if(pos>0)
	    {
            pixel.push_back(count);
	    }
        count++;
    }

    Eigen::MatrixXd Matrix_M(n_abs_4, n_abs_4);
    Matrix_M = Eigen::MatrixXd::Zero(n_abs_4, n_abs_4);
	Matrix_M(0, 0) = 1 / t_el;
	Matrix_M(0, 1) = L_I * InputPara.G_tes_bath / (I * InputPara.L);

	Matrix_M(1, 0) = -I * InputPara.R * (2 + InputPara.beta) / InputPara.C_tes;
	Matrix_M(1, 1) = 1 / t_I + (InputPara.G_abs_tes / InputPara.C_tes);
	Matrix_M(1, 2) = -InputPara.G_abs_tes / InputPara.C_tes;

	Matrix_M(2, 1) = -InputPara.G_abs_tes/C_abs;
	Matrix_M(2, 2) = InputPara.G_abs_tes / C_abs+G_abs_abs/C_abs;
	Matrix_M(2, 3) = -G_abs_abs/C_abs;

	Matrix_M(n_abs_1, InputPara.n_abs) = -G_abs_abs / C_abs;
	Matrix_M(n_abs_1, n_abs_1) = InputPara.G_abs_tes / C_abs + G_abs_abs / C_abs;
	Matrix_M(n_abs_1, n_abs_2) = -InputPara.G_abs_tes / C_abs;

	Matrix_M(n_abs_2, n_abs_1) = -InputPara.G_abs_tes / InputPara.C_tes;
	Matrix_M(n_abs_2, n_abs_2) = 1 / t_I + InputPara.G_abs_tes / InputPara.C_tes;
	Matrix_M(n_abs_2, n_abs_3) = -I * InputPara.R * (2 + InputPara.beta) / InputPara.C_tes;

	Matrix_M(n_abs_3, n_abs_2) = L_I * InputPara.G_tes_bath / (I * InputPara.L);
    Matrix_M(n_abs_3, n_abs_3) = 1 / t_el;
	for(int i=0;i<=n_abs_4;i++)
	{
		for(int j;j<=n_abs_4;j++)
		{
			if(j==i-1)
			{
				Matrix_M(i, j) = -G_abs_abs / C_abs;
			}
			if(j==i)
			{
				Matrix_M(i, j) = 2 * G_abs_abs / C_abs;
			}
			if(j==i+1)
			{
				Matrix_M(i, j) = -G_abs_abs / C_abs;
			}
		}
	}
    Matrix_M *= -1;

    Eigen::MatrixXd Matrix_X(n_abs_2,n_abs_4);
    Matrix_X = Eigen::MatrixXd::Zero(n_abs_2, n_abs_4);

    const double Same = InputPara.E * 1e3 * e / InputPara.C_tes;
    Matrix_X(0, 1) = Same;
    Matrix_X(n_abs_1,n_abs_2)= Same;
    for(int i=0;i<=n_abs_4;i++)
    {
        for (const int& pix : pixel)
        {
            if(i==pix)
            {
                Matrix_X(i,i+1) = InputPara.positions[i-1];
            }
        }
    }
    Eigen::EigenSolver<Eigen::MatrixXd> solver(Matrix_M);
    // 固有値
    Eigen::VectorXcd EigenValues = solver.eigenvalues();
    // 固有ベクトル
    Eigen::MatrixXcd EigenVectors = solver.eigenvectors();
    std::vector<Eigen::VectorXcd> arb;
    
    // Xの各列に対して線形方程式を解く
    for (int i = 0; i < Matrix_X.cols(); ++i) {
        Eigen::VectorXcd right_hand_side = Matrix_X.col(i); // Xのi番目の列を右辺に設定
        Eigen::VectorXcd solution = EigenVectors.lu().solve(right_hand_side); // 線形方程式を解く
        arb.push_back(solution); // 結果をarbに追加
    }
    std::vector<double> time= linspace(0,InputPara.samples/InputPara.rate, static_cast<int>(InputPara.samples));

    double I_t0_sum = 0;
    double I_t1_sum = 0;
    count = 0;
    I_t0.resize(pixel.size(), time.size());
    std::cout << "pixel:" << pixel.size() << "\ntime:" << time.size() << "\n";
    for(const int& pix:pixel)
    {
        Eigen::MatrixXcd Matrix_t(n_abs_4,time.size());
        Matrix_t.setZero();
        std::cout << "Model 4.1\n";
	    for(int j=0;j<=n_abs_4;j++)
	    {
            int counter = 0;
            //std::cout << "result:" << arb[pix](j) * EigenVectors(0, j) * std::exp(EigenValues(j)*time[0]) << "\n";
            std::cout << "result:" << arb[pix](j);
            std::cout <<"Matrix:"<< Matrix_t(j, 0) << "\n";
            std::cout << "counter:" << counter << "\n";
            for(const auto& ti:time)
            {
               //std::cout << "j:" << j << "\ncounter" << counter << "\n";
                Matrix_t(j,counter) = arb[pix](j) * EigenVectors(0, j) * std::exp(EigenValues(j) * ti);
                counter++;
            }
            
	    }
        std::cout << "Model 4.4\n";
        I_t0.row(count) = Matrix_t.colwise().sum().real();
        count++;
    }
    std::cout << "Model5\n";
    count = 0;
    I_t1.resize(pixel.size(), time.size());
    for (const int& pix : pixel)
    {
        Eigen::MatrixXcd Matrix_t;
        for (int j = 0; j <= n_abs_4; j++)
        {
            int counter = 0;
            for (const auto& ti : time)
            {
                Matrix_t(j, counter) = arb[pix](j) * EigenVectors(n_abs_3, j) * std::exp(EigenValues(j) * ti);
                counter++;
            }
        }
        Eigen::VectorXd ColSum = Matrix_t.colwise().sum().real();
        count++;
        I_t1.row(count)=ColSum;
    }

    I_t0 *= -1;
    I_t1 *= -1;
    std::cout << "Model 6\n";
}
