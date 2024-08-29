#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <filesystem>
#include <sstream>
#include <unordered_set>
#include <nlohmann/json.hpp>

#include"Dump2Batch.hpp"

//#define Python



#ifdef Python
extern "C" __declspec(dllexport) void MakeOutput(const char* DataPath, const char* InputPath_char, const char* Output_FileName) {

    std::string path(DataPath);
    std::string InputPath(InputPath_char);
    std::string output_file(Output_FileName);
#else
void main(){
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
        
    std::map<int, std::map<int, EventInfo>> batch = ReadDump(DumpPath);

    std::cout<<"Finished\nWriting output.json...\n";

    WriteOutput(batch, output_file);

    std::cout << "Finished\n";

    return;
}
