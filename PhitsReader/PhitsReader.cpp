#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <nlohmann/json.hpp>
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
    std::string DataPath = "./";
    std::string InputPath = "./input.json";
    std::string output_file = "output.json";
#endif
    InputParameters InputPara=ReadInputJson(InputPath);
    std::cout << "Input.json is parsed\n";
#ifdef Python
    DataPath += ("/" + InputPara.output+"/");
#endif


    std::string DumpPath = DataPath + "dumpall.dat";

    std::cout << "Processing file...\n";
        
    std::map<int, std::map<int, EventInfo>> batch;
	ReadDump(DumpPath, batch);

    std::cout<<"Finished\nWriting output.json...\n";

    WriteOutput(batch, output_file);

    PulseParameters PulsePara(InputPara);

    std::cout << "Finished\n";

    Eigen::MatrixXd Matrix_M=MakeMatrix_M(PulsePara, InputPara);

    Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(Matrix_M);
    // 固有値
    Eigen::VectorXcd EigenValues = eigensolver.eigenvalues();
    // 固有ベクトル
    Eigen::MatrixXcd EigenVectors = eigensolver.eigenvectors();

    // ファイルへの書き出し
    std::ofstream file("eigenvalues.txt");
    if (file.is_open()) {
        for (int i = 0; i < EigenValues.size(); ++i) {
            file << EigenValues(i).real() << "\n";
        }
        file.close();
    }
    else {
        std::cerr << "Unable to open file for writing.\n";
    }

    return 0;
}
