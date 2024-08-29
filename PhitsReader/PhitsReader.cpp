#pragma once
#include <iostream>

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

    return 0;
}
