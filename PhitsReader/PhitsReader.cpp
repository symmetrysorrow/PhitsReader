﻿#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <nlohmann/json.hpp>
#include <iomanip>
#include <thread>
#include <ppl.h> 
#include "SpinProgress.h"
#include "Batch2Pulse.h"
#include"Dump2Batch.h"

//#define Python

bool DEBUG = true;

#ifdef Python
extern "C" __declspec(dllexport) int MakeOutput(const char* DataPath_char, const char* InputPath_char) {

    std::string DataPath(DataPath_char);
    std::string InputPath(InputPath_char);
    //std::string output_file(Output_FileName);
#else
int main(){
    std::string DataPath = ".";
    std::string InputPath = "./input.json";
    std::string output_file = "output.json";
#endif
    InputParameters InputPara=ReadInputJson(InputPath);
    std::cout << "Input.json is parsed\n";
#ifdef Python
    DataPath += ("/" + InputPara.output);
#endif


    std::string DumpPath = DataPath + "/dumpall.dat";

	SpinProgress spinner_braille(SpinProgress::Braille);
	spinner_braille.set_message("Processing dump file...");
        
    std::map<int, std::map<int, EventInfo>> batch;
	int ReadReturn = ReadDump(DumpPath, batch);

	if(ReadReturn==-1)
	{
		spinner_braille.complete("Error!");
		std::cerr << "Failed to open dump file: " << DumpPath << std::endl;
		return -1;
	}

	spinner_braille.complete("Completed");

    //WriteOutput(batch, output_file);

    PulseParameters PulsePara(InputPara);
	
	// 出力ディレクトリを作成
	std::filesystem::create_directories(DataPath + "/Pulse/Ch0");
	std::filesystem::create_directories(DataPath + "/Pulse/CH1");

    Eigen::MatrixXd Matrix_M=MakeMatrix_M(PulsePara, InputPara);

	Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(Matrix_M);
	// 固有値
	Eigen::VectorXcd EigenValues = eigensolver.eigenvalues();
	// 固有ベクトル
	Eigen::MatrixXcd EigenVectors = eigensolver.eigenvectors();
	
	//SortEigen(EigenValues, EigenVectors);

	const int n_abs = InputPara.n_abs;
	const int n_abs_1 = n_abs + 1;
	const int n_abs_2 = n_abs + 2;
	const int n_abs_3 = n_abs + 3;
	const int n_abs_4 = n_abs + 4;

	const std::vector<double> Block = linspace(-1, 1, n_abs + 1);

	// 全体のサイズ
	int total_size = static_cast<int>(batch.size());
	// 処理した要素の数をカウントする原子変数
	std::atomic<int> processed_count(0);
	// 表示した進捗パーセンテージを追跡する
	std::atomic<int> last_percent_displayed(0);

	concurrency::parallel_for_each(batch.begin(), batch.end(), [&](const std::pair<const int, std::map<int, EventInfo>>& outer_pair) {
		
		std::vector<double> BlockDeposit(n_abs, 0.0);

		std::vector<double> Positions;

		for (const auto& inner_pair : outer_pair.second) {
			for (int i = 0; i < inner_pair.second.x_deposit.size(); i++)
			{
				int Pixel = InBlock(Block, inner_pair.second.x_deposit[i], inner_pair.second.y_deposit[i], inner_pair.second.z_deposit[i]);
				BlockDeposit[Pixel - 1] += inner_pair.second.E_deposit[i];
			}
		}
		for (auto& block_deposit : BlockDeposit)
		{
			Positions.push_back(block_deposit * 1e6);
		}

		std::vector<int> pixel;
		int count = 1;

		for (const double& pos : Positions)
		{
			if (pos > 0)
			{
				pixel.push_back(count);
			}
			count++;
		}

		Eigen::MatrixXd Matrix_X = MakeMatrix_X(PulsePara, InputPara, pixel, Positions);
		

		Eigen::MatrixXd EiVec = EigenVectors.real();

		// arb行列の初期化
		Eigen::MatrixXd arb(Matrix_X.rows(), Matrix_X.cols());

		arb.setZero();

		// 各ベクトルごとに計算を行い、結果を arb に格納します。
		for (int i = 0; i < Matrix_X.rows(); ++i) {
			Eigen::VectorXd b = Matrix_X.row(i); // 行ベクトルを取得
			Eigen::VectorXd x = EiVec.colPivHouseholderQr().solve(b); // 連立方程式を解く
			arb.row(i) = x;

		}

		std::vector<double> time = linspace(0, InputPara.samples / InputPara.rate, static_cast<int>(InputPara.samples));

		Eigen::MatrixXd pulse_total_0;
		Eigen::MatrixXd pulse_total_1;

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
		}
		pulse_total_0 *= -1;

		pulse_total_1.resize(pixel.size(), time.size());

		for (int i = 0; i < pixel.size(); i++)
		{
			Eigen::MatrixXcd Matrix_t(n_abs_4, time.size());
			Matrix_t.setZero();
			for (int j = 0; j < n_abs_4; j++)
			{
				for (int k = 0; k < time.size(); k++)
				{
					Matrix_t(j, k) = arb(pixel[i], j) * EigenVectors(n_abs_3, j) * std::exp(EigenValues(j) * time[k]);
				}
			}
			pulse_total_1.row(i) = Matrix_t.colwise().sum().real();
		}

		pulse_total_1 *= -1;

		Eigen::VectorXd pulse_0 = pulse_total_0.colwise().sum();
		Eigen::VectorXd pulse_1 = pulse_total_1.colwise().sum();

		std::string ChFile_0 = DataPath + "/Pulse/Ch0/CH0_" + std::to_string(outer_pair.first) + ".dat";
		std::ofstream outFile_0(ChFile_0);
		if (!outFile_0) {
			std::cerr << "Failed to open file:" <<ChFile_0<< std::endl;
			return -1;
		}

		for (int i = 0; i < pulse_0.size(); ++i) {
			outFile_0 << pulse_0(i) << std::endl;
		}

		outFile_0.close();

		std::string ChFile_1 = DataPath + "/Pulse/Ch1/CH1_" + std::to_string(outer_pair.first) + ".dat";
		std::ofstream outFile_1(ChFile_1);
		if (!outFile_1) {
			std::cerr << "Failed to open file:" <<ChFile_1<< std::endl;
			return -1;
		}

		for (int i = 0; i < pulse_1.size(); ++i) {
			outFile_1 << pulse_1(i) << std::endl;
		}

		outFile_1.close();

		// 処理した要素のカウントをインクリメント
		int current_count = ++processed_count;
		// 進捗パーセンテージの計算
		int current_percent = (current_count * 100) / total_size;
		// 進捗パーセンテージの表示（1%単位で更新）
		if (current_percent > last_percent_displayed) {
			last_percent_displayed = current_percent;
			std::cout << current_percent << "% Converting to pulse...\r" << std::flush;
		}
	});
	// 進捗表示行を消去して、結果を表示
	std::cout << "\r\033[K";  // カーソルを行の先頭に戻し、行をクリア
	std::cout << u8"✓ Completed" << std::endl;
	std::cout << "Finished!";
    return 0;
}
