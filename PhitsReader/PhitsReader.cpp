#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <nlohmann/json.hpp>
#include <iomanip>
#include <thread>
#include <ppl.h> 
#include <mutex>
#include "Batch2Pulse.h"
#include"Dump2Batch.h"
#include <concurrent_vector.h>
#include"pulse2csv.h"
#include "AddNoise.h"

int main(){

	std::string DataPath;
	std::cout << "Data Path:";
	std::cin >> DataPath;
    std::string InputPath = DataPath+"/input.json";

    InputParameters InputPara=ReadInputJson(InputPath);

	std::vector<std::string> DumpPathes;

	

	for (const auto& posi : InputPara.positions) {
		DumpPathes.push_back(DataPath+"\\"+ std::to_string(static_cast<int>(InputPara.output * 1000))+"keV_"+std::to_string(posi));
	}
	int TotalCounter = 1;
	for (const auto& DumpPath : DumpPathes) {
		std::cout << TotalCounter << "/" << DumpPathes.size() <<"\n";

		std::string DumpFilePath = DumpPath + "\\dumpall.dat";

		std::cout << "Processing file...\n";

		std::map<int, std::map<int, EventInfo>> batch;
		int ReadReturn = ReadDump(DumpFilePath, batch);

		if (ReadReturn == -1)
		{
			return -1;
		}

		std::cout << "Finished\n";

		std::cout << "Converting to pulse...\n";

		PulseParameters PulsePara(InputPara);
#if 0
		// 出力ディレクトリを作成
		std::filesystem::create_directories(DataPath + "/Pulse/Ch0");
		std::filesystem::create_directories(DataPath + "/Pulse/CH1");
#endif
		Eigen::MatrixXd Matrix_M = MakeMatrix_M(PulsePara, InputPara);

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

		int Counter = 0;

		concurrency::concurrent_vector<std::tuple<int, double, double>> PulseInfo_Ch0;
		concurrency::concurrent_vector<std::tuple<int, double, double>> PulseInfo_Ch1;

		std::pair<std::vector<double>, std::vector<double>> Coeffs = MakeCoeff(InputPara);
		Eigen::VectorXd Noise_dense = readLinesToEigen("./noise_spectral_total_alpha71beta1.6.dat");

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
					for (int k = 0; k < time.size(); k++)
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

			if (InputPara.noise) {
				AddNoise(Noise_dense, pulse_0);
				AddNoise(Noise_dense, pulse_1);
			}

			PulseInfo_Ch0.push_back(MakeCSV(pulse_0, outer_pair.first, Coeffs.first, Coeffs.second));
			PulseInfo_Ch1.push_back(MakeCSV(pulse_1, outer_pair.first, Coeffs.first, Coeffs.second));

		});

		std::string ChFile_0 = DumpPath + "/output_TES0_with_noise.csv";
		std::ofstream outFile_0(ChFile_0);
		if (!outFile_0) {
			std::cerr << "Failed to open file:" << ChFile_0 << std::endl;
			return -1;
		}
		outFile_0 << ", height,rise\n";
		for (const std::tuple<int, double, double>& info : PulseInfo_Ch0) {
			outFile_0 << std::get<0>(info) << "," << std::get<1>(info) << "," << std::get<2>(info) << "\n";
		}
		outFile_0.close();

		std::string ChFile_1 = DumpPath + "/output_TES1_with_noise.csv";
		std::ofstream outFile_1(ChFile_1);
		if (!outFile_1) {
			std::cerr << "Failed to open file:" << ChFile_1 << std::endl;
			return -1;
		}
		outFile_1 << ",height,rise\n";
		for (const std::tuple<int, double, double>& info : PulseInfo_Ch1) {
			outFile_1 << std::get<0>(info) << "," << std::get<1>(info) << "," << std::get<2>(info) << "\n";
		}
		outFile_1.close();
		TotalCounter++;
	}

	std::cout << "Finished!\n";
    
    return 0;
}
