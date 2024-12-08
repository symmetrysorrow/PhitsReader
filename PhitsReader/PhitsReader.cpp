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
#include "SpinProgress.h"

int main(){

	std::string DataPath;
	std::cout << "Data Path:";
	std::cin >> DataPath;
    std::string InputPath = DataPath+"/input.json";

    InputParameters InputPara=ReadInputJson(InputPath);

	std::vector<std::string> DumpPathes;

	for (const auto& posi : InputPara.positions) {
		DumpPathes.push_back(DataPath+"/"+ std::to_string(static_cast<int>(InputPara.E))+"keV_"+std::to_string(posi));
	}

	int TotalCounter = 1;

	for (const auto& DumpPath : DumpPathes) {
		std::cout << TotalCounter << "/" << DumpPathes.size() << "\n";

		SpinProgress spinner;
		spinner.set_message("Processing dumpall file...");

		std::string DumpFilePath = DumpPath + "/dumpall.dat";

		std::map<int, std::map<int, EventInfo>> batch;
		int ReadReturn = ReadDump(DumpFilePath, batch,InputPara.E/1000);

		if (ReadReturn == -1)
		{
			std::cout << "Error in dump file\n";
			return -1;
		}

		spinner.complete("Finished");

		PulseParameters PulsePara(InputPara);
		std::filesystem::create_directories(DumpPath + "/Pulse_ms/Ch0");
		std::filesystem::create_directories(DumpPath + "/Pulse_ms/CH1");

		std::filesystem::create_directories(DumpPath + "/Pulse_ms_noise/Ch0");
		std::filesystem::create_directories(DumpPath + "/Pulse_ms_noise/CH1");

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

		const std::vector<double> Block = linspace(-InputPara.length/20, InputPara.length/20, n_abs + 1);

		int Counter = 0;

		concurrency::concurrent_vector<std::tuple<int, double, double>> PulseInfo_Ch0;
		concurrency::concurrent_vector<std::tuple<int, double, double>> PulseInfo_Ch1;

		std::pair<std::vector<double>, std::vector<double>> Coeffs = MakeCoeff(InputPara);
		std::string NoisePath = DataPath + "/noise_spectral_total_alpha71beta1.6.dat";
		Eigen::VectorXd Noise_dense = readLinesToEigen(NoisePath);

		size_t total_items = batch.size(); // 全体の要素数を取得
		std::atomic<size_t> completed_items(0); // 完了したアイテム数
		std::mutex output_mutex; // 出力用ミューテックス

		concurrency::parallel_for_each(batch.begin(), batch.end(), [&](const std::pair<const int, std::map<int, EventInfo>>& outer_pair) {

			std::vector<double> BlockDeposit(n_abs, 0.0);

			std::vector<double> Positions;

			for (const auto& inner_pair : outer_pair.second) {
				for (int i = 0; i < inner_pair.second.x_deposit.size(); i++)
				{
					int Pixel = InBlock(InputPara,Block, inner_pair.second.x_deposit[i], inner_pair.second.y_deposit[i], inner_pair.second.z_deposit[i]);
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

			std::string PulseFile_0 = DumpPath + "/Pulse_ms/Ch0/CH0_" + std::to_string(outer_pair.first) + ".dat";
			std::string PulseFile_1 = DumpPath + "/Pulse_ms/Ch1/CH1_" + std::to_string(outer_pair.first) + ".dat";

			std::ofstream PulseoutFile_0(PulseFile_0);
			if (!PulseoutFile_0) {
				std::cerr << "Failed to open file:" << PulseFile_0 << std::endl;
				return -1;
			}
			for (int i = 0; i < pulse_0.size(); ++i) {
				PulseoutFile_0 << pulse_0[i] << std::endl;
			}
			PulseoutFile_0.close();

			std::ofstream PulseoutFile_1(PulseFile_1);
			if (!PulseoutFile_1) {
				std::cerr << "Failed to open file:" << PulseFile_1 << std::endl;
				return -1;
			}
			for (int i = 0; i < pulse_0.size(); ++i) {
				PulseoutFile_1 << pulse_1[i] << std::endl;
			}
			PulseoutFile_1.close();


			AddNoise(Noise_dense, pulse_0,InputPara.rate);
			AddNoise(Noise_dense, pulse_1, InputPara.rate);

			PulseFile_0 = DumpPath + "/Pulse_ms_noise/Ch0/CH0_" + std::to_string(outer_pair.first) + ".dat";
			PulseFile_1 = DumpPath + "/Pulse_ms_noise/Ch1/CH1_" + std::to_string(outer_pair.first) + ".dat";

			std::vector<double> Filterd_CH0 = ApplyFilter(pulse_0, Coeffs.first, Coeffs.second);
			std::vector<double> Filterd_CH1 = ApplyFilter(pulse_1, Coeffs.first, Coeffs.second);

			PulseInfo_Ch0.push_back(GetPulseInfo(outer_pair.first,Filterd_CH0));
			PulseInfo_Ch1.push_back(GetPulseInfo(outer_pair.first,Filterd_CH1));

			std::ofstream PulseoutFile_0_n(PulseFile_0);
			if (!PulseoutFile_0_n) {
				std::cerr << "Failed to open file:" << PulseFile_0 << std::endl;
				return -1;
			}
			for (int i = 0; i < Filterd_CH0.size(); ++i) {
				PulseoutFile_0_n << Filterd_CH0[i] << std::endl;
			}
			PulseoutFile_0_n.close();

			std::ofstream PulseoutFile_1_n(PulseFile_1);
			if (!PulseoutFile_1_n) {
				std::cerr << "Failed to open file:" << PulseFile_1 << std::endl;
				return -1;
			}
			for (int i = 0; i < Filterd_CH1.size(); ++i) {
				PulseoutFile_1_n << Filterd_CH1[i] << std::endl;
			}
			PulseoutFile_1_n.close();

			size_t completed = completed_items.fetch_add(1);

			// 進捗を更新（1%ごと）
			size_t progress = static_cast<size_t>((completed + 1) * 100 / total_items); // +1は現在のアイテムを含むため
			if (progress > 0 && progress % 1 == 0) { // 1%ごとに更新
				std::lock_guard<std::mutex> lock(output_mutex); // スレッドセーフな出力
				std::cout << "\rConverting to Pulse: " << progress << "%" << std::flush; // プログレスを表示
			}
		});

		{
			std::lock_guard<std::mutex> lock(output_mutex);
			std::cout << "\nFinished\n";
		}
		TotalCounter++;
	}

	std::cout << "Completed!\n";
    
    return 0;
}
