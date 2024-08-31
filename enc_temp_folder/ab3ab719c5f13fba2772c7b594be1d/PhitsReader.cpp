#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <nlohmann/json.hpp>
#include <iomanip>
#include "Batch2Pulse.h"
#include"Dump2Batch.h"

//#define Python

#ifdef Python
extern "C" __declspec(dllexport) void MakeOutput(const char* DataPath, const char* InputPath_char, const char* Output_FileName) {

    std::string path(DataPath);
    std::string InputPath(InputPath_char);
    std::string output_file(Output_FileName);
#else
int main(){
    std::string DataPath = ".";
    std::string InputPath = "./input.json";
    std::string output_file = "output.json";
#endif
    InputParameters InputPara=ReadInputJson(InputPath);
    std::cout << "Input.json is parsed\n";
#ifdef Python
    DataPath += ("/" + InputPara.output+"/");
#endif


    std::string DumpPath = DataPath + "/dumpall.dat";

    std::cout << "Processing file...\n";
        
    std::map<int, std::map<int, EventInfo>> batch;
	ReadDump(DumpPath, batch);

    std::cout<<"Finished\nWriting output.json...\n";

    WriteOutput(batch, output_file);

    PulseParameters PulsePara(InputPara);

    std::cout << "Finished\n";

	// 出力ディレクトリを作成
	std::filesystem::create_directories(DataPath + "/PulseCPP/Ch0");
	std::filesystem::create_directories(DataPath + "/PulseCPP/CH1");

    Eigen::MatrixXd Matrix_M=MakeMatrix_M(PulsePara, InputPara);

	Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(Matrix_M);
	// 固有値
	Eigen::VectorXcd EigenValues = eigensolver.eigenvalues();
	// 固有ベクトル
	Eigen::MatrixXcd EigenVectors = eigensolver.eigenvectors();

	SortEigen(EigenValues, EigenVectors);

	std::ofstream file_m("vec.txt");

	if (file_m.is_open()) {

		file_m << std::scientific << std::setprecision(25) << EigenVectors.real() << std::endl;

		file_m.close();
	}

	std::ofstream file_val("val.txt");

	if (file_val.is_open()) {

		file_val << std::scientific << std::setprecision(25) << EigenValues.real() << std::endl;

		file_val.close();
	}

#if 0
	
	// 固有値と固有ベクトルが元の行列に対して正しいか検算する
	for (int i = 0; i < EigenVectors.cols(); ++i) {
		Eigen::VectorXd v = EigenVectors.col(i).real();
		double lambda = EigenValues(i).real();

		// A * v を計算
		Eigen::VectorXd Av = Matrix_M * v;

		// lambda * v を計算
		Eigen::VectorXd lambda_v = lambda * v;

		// 比較結果を表示
		if ((Av - lambda_v).norm() < 1e-5) {  // 許容誤差を設定
			std::cout << "Eigenvalue and Eigenvector are correct for index " << i << "." << std::endl;
		}
		else {
			std::cout << "Eigenvalue and Eigenvector are incorrect for index " << i << "." << std::endl;
			std::cout << "A * v: " << Av.transpose() << std::endl;
			std::cout << "lambda * v: " << lambda_v.transpose() << std::endl;
		}
	}
#endif


	const int n_abs = InputPara.n_abs;
	const int n_abs_1 = n_abs + 1;
	const int n_abs_2 = n_abs + 2;
	const int n_abs_3 = n_abs + 3;
	const int n_abs_4 = n_abs + 4;

	const std::vector<double> Block = linspace(-1, 1, n_abs + 1);

	int Counter = 0;
	int sum = batch.size();
	for (const auto& outer_pair : batch) {
		{
			float progress = Counter == sum ? 100 : static_cast<float>(Counter) / sum * 100;
			std::cout << "\rProgress:" << std::fixed << std::setprecision(2) << progress << "%";
			std::cout.flush();
			Counter++;
		}

		std::vector<double> BlockDeposit(n_abs, 0.0);

		InputPara.positions.clear();

		for (const auto& inner_pair : outer_pair.second) {
			for (int i = 0; i < inner_pair.second.x_deposit.size(); i++)
			{
				int Pixel = InBlock(Block, inner_pair.second.x_deposit[i], inner_pair.second.y_deposit[i], inner_pair.second.z_deposit[i]);
				BlockDeposit[Pixel - 1] += inner_pair.second.E_deposit[i];
				//std::cout << "E_dep[" << i << "]:" << std::scientific << std::setprecision(15) << inner_pair.second.E_deposit[i] << "\n";

			}
		}
		for (auto& block_deposit : BlockDeposit)
		{
			InputPara.positions.push_back(block_deposit * 1e6);
		}

		std::vector<int> pixel;
		int count = 1;

		for (const double& pos : InputPara.positions)
		{
			if (pos > 0)
			{
				pixel.push_back(count);
			}
			count++;
		}

		Eigen::MatrixXd Matrix_X = MakeMatrix_X(PulsePara, InputPara, pixel);

		Eigen::MatrixXd EiVec = EigenVectors.real();

		// arb行列の初期化
		Eigen::MatrixXd arb(Matrix_X.rows(), Matrix_X.cols());

		// 各ベクトルごとに計算を行い、結果を arb に格納します。
		for (int i = 0; i < Matrix_X.rows(); ++i) {
			Eigen::VectorXd b = Matrix_X.row(i); // 行ベクトルを取得
			Eigen::VectorXd x = EiVec.colPivHouseholderQr().solve(b); // 連立方程式を解く

			arb.row(i) = x;
		}

		std::vector<double> time = linspace(0, InputPara.samples / InputPara.rate, static_cast<int>(InputPara.samples));

		Eigen::MatrixXd pulse_total_0;
		Eigen::MatrixXd pulse_total_1;

		count = 0;
		pulse_total_0.resize(pixel.size(), time.size());

		for (int i = 0; i < pixel.size(); i++)
		{
			Eigen::MatrixXcd Matrix_t(n_abs_4, time.size());
			Matrix_t.setZero();
			for (int j = 0; j < n_abs_4; j++)
			{
				for (int k=0;k<time.size();k++)
				{
					Matrix_t(j, k) = arb(pixel[i], j) * EigenVectors(0, j) * std::exp(EigenValues(j) * time[k]);
				}
			}
			pulse_total_0.row(i) = Matrix_t.colwise().sum().real();

			std::ofstream file_t("t.txt");
			if (file_t.is_open()) {

				file_t << Matrix_t.real() << std::endl;

				file_t.close();
				int something = 0;
				std::cin >> something;
				something += 1;
			}
			count++;
		}
#if 0
		std::ofstream file("arb.txt");

		if (file.is_open()) {

			file << arb << std::endl;

			file.close();
		}
		std::ofstream file_x("x.txt");
		if (file_x.is_open()) {

			file_x << Matrix_X << std::endl;

			file_x.close();
		}
#endif
		count = 0;
		pulse_total_1.resize(pixel.size(), time.size());

		for (int i=0;i<pixel.size();i++)
		{
			Eigen::MatrixXcd Matrix_t(n_abs_4, time.size());
			Matrix_t.setZero();
			for (int j = 0; j < n_abs_4; j++)
			{
				int counter = 0;
				for (const auto& ti : time)
				{
					Matrix_t(j, counter) = arb(pixel[i], j) * EigenVectors(n_abs_3, j) * std::exp(EigenValues(j) * ti);
					counter++;
				}
			}
			pulse_total_1.row(count) = Matrix_t.colwise().sum().real();
			
			count++;
		}

		pulse_total_0 *= -1;
		pulse_total_1 *= -1;

		Eigen::VectorXd pulse_0 = pulse_total_0.colwise().sum();
		Eigen::VectorXd pulse_1 = pulse_total_1.colwise().sum();

		std::string ChFile_0 = DataPath + "/PulseCPP/Ch0/CH0_" + std::to_string(outer_pair.first) + ".dat";
		std::ofstream outFile_0(ChFile_0);
		if (!outFile_0) {
			std::cerr << "ファイルを開けませんでした。" << std::endl;
			return -1;
		}

		for (int i = 0; i < pulse_0.size(); ++i) {
			outFile_0 << pulse_0(i) << std::endl;
		}

		outFile_0.close();

		std::string ChFile_1 = DataPath + "/PulseCPP/Ch1/CH1_" + std::to_string(outer_pair.first) + ".dat";
		std::ofstream outFile_1(ChFile_1);
		if (!outFile_1) {
			std::cerr << "ファイルを開けませんでした。" << std::endl;
			return -1;
		}

		for (int i = 0; i < pulse_1.size(); ++i) {
			outFile_1 << pulse_1(i) << std::endl;
		}

		outFile_1.close();
	}
	
    return 0;
}
