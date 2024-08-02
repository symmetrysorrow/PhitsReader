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
#include"SubFuncs.hpp"
#include <Eigen/Dense>
//#define Python



#ifdef Python
extern "C" __declspec(dllexport) void MakeOutput(const char* DataPath_char, const char* InputPath_char) {

    std::string DataPath(DataPath_char);
    std::string InputPath(InputPath_char);
#else
void main(){
    std::string DataPath = "F:/hata/output_5";
    std::string InputPath="./input.json";
#endif

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

    InputParameters InputPara;

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

    DataPath += ("/"+InputPara.output);
    
    //定数パラメーター
    constexpr double emin_electron = 0.1;
    constexpr double emin_photon = 0.001;

    //all plot[0], one plot[eventnumber], no plot[-1]
    constexpr float event_number = -1;

    //各イベントの情報をhistoryに代入し、適宜batchに入力する。最終結果はbatchに入る。
    std::map<int, EventInfo> history;
    std::map<int, std::map<int, EventInfo>> batch;

    //計算に使われる変数
    double ncol = 1;
    std::vector<double> xyz= {0,0,0};
    std::vector<double> cxyz = {0,0,0};
    double cnt = 0;
    double num = 0;
    double nocas = 0;
    double no = 0;
    std::vector<double> reg(2);
    std::vector<double> name;
    double benergy, cenergy, ityp, nclsts, jcoll, energy_new, ncl, energy, energy_dps;

    std::ifstream DumpFile(DataPath+"/dumpall.dat", std::ios::binary);

    if (!DumpFile.is_open()) {
        std::cerr << "Failed to open file: " << DataPath << std::endl;
        return;
    }

    std::cout << "Processing file...\n";

    std::string line;
    while (std::getline(DumpFile, line)) //ファイルの各行ごとに実行
    {
        std::vector<double> column=split_line(line);
        
        if (static_cast<int>(ncol) == 1)
        {
            ncol = 4;
            cnt = 0;
            nocas = 0;
            no = 0;
            continue;
        }
        std::unordered_set<int> valid_values = { 1, 2, 3, 17 };
        if (valid_values.find(static_cast<int>(ncol)) == valid_values.end())
        {
            if (static_cast<int>(cnt) == 0)
            {
                ncol = column[0];
                if (ncol != 4) { cnt += 1; }
                else
                {
                    if (no > 1)
                    {
                        batch[static_cast<int>(nocas)] = history;
                        history.clear();
                        std::map<int, EventInfo> emptyMap;
                        history.swap(emptyMap);

                        if (event_number == nocas || static_cast<int>(event_number) == 0) {}
                    }
                    else {}

                    EventInfo event;
                    event.ityp = 14;
                    history[1] = event;
                }
            }

            if (static_cast<int>(cnt) == 1 && static_cast<int>(ncol) == 4) {nocas = column[0]; }

            if (static_cast<int>(cnt) == 2)
            {
                no = column[0];
                ityp = column[2];
                if (ityp != 12 && ityp != 13) { cnt += 1; }
            }

            if (static_cast<int>(cnt) == 5 && reg.size()>=2) {std::copy_n(column.begin(), 2, reg.begin()); }

            if (static_cast<int>(cnt) == 8) { name = column; }

            if (static_cast<int>(cnt) == 11)
            {
                benergy = column[0];
                xyz[1] = column[0];
                xyz[2] = column[1];
            }

            if (static_cast<int>(cnt) == 13)
            {
                cenergy = column[0];
                cxyz[0] = column[2];
            }

            if (static_cast<int>(cnt) == 14)
            {
                cxyz[1] = column[0];
                cxyz[2] = column[1];
            }
            if (static_cast<int>(cnt) == 16)
            {
                if (!(static_cast<int>(ncol) == 13 || static_cast<int>(ncol) == 14)) { cnt = -1; }
                if (static_cast<int>(ityp) == 14 || static_cast<int>(ityp) == 12 || static_cast<int>(ityp) == 13) {
                    if (ncol == 4) { }
                   
                    if (history.find(static_cast<int>(no)) == history.end())
                    {
                        EventInfo new_event;
                        new_event.ityp = static_cast<int>(ityp);
                        history[static_cast<int>(no)] = new_event;
                    }
                    history[static_cast<int>(no)].x.push_back(cxyz[0]);
                    history[static_cast<int>(no)].y.push_back(cxyz[1]);
                    history[static_cast<int>(no)].z.push_back(cxyz[2]);
                    history[static_cast<int>(no)].E.push_back(cenergy);

                    if (static_cast<int>(ncol) == 11)
                    {
                        history[static_cast<int>(no)].E_deposit.push_back(energy);
                        history[static_cast<int>(no)].x_deposit.push_back(cxyz[0]);
                        history[static_cast<int>(no)].y_deposit.push_back(cxyz[1]);
                        history[static_cast<int>(no)].z_deposit.push_back(cxyz[2]);
                        
                    }
                }
            }

            if (static_cast<int>(cnt) == 17) { nclsts = column[0]; }

            if (static_cast<int>(cnt) == 18)
            {
                jcoll = column[2];
                ncol = 17;
                cnt = -1;
                ncl = 0;
                energy_new = 0;

                if (static_cast<int>(jcoll) == 14)
                {
                    cnt += 1;
                    num += 1;
                    continue;
                }
            }

        }

        if (static_cast<int>(ncol) == 17)
        {
            if (static_cast<int>(cnt == 1)) { ityp = column[3]; }
            if (static_cast<int>(cnt) == 5) { energy = column[1]; }
            if (static_cast<int>(cnt) == 8)
            {
	            if (int iityp = static_cast<int>(ityp); iityp == 14 || iityp == 12 || iityp == 13)
                {
                    if ((energy < emin_electron && iityp == 12) || (energy >= emin_photon && iityp == 14) || (energy >= emin_electron && iityp == 13))
                    {
                        energy_new += energy;
                    }
                    energy_dps = benergy - energy_new;

                    if (static_cast<int>(ncl) == static_cast<int>(nclsts) - 1)
                    {
                        ncol = 13;
                        history[static_cast<int>(no)].E_deposit.push_back(energy_dps);
                        history[static_cast<int>(no)].x_deposit.push_back(cxyz[0]);
                        history[static_cast<int>(no)].y_deposit.push_back(cxyz[1]);
                        history[static_cast<int>(no)].z_deposit.push_back(cxyz[2]);
                    }
                    else { ncl += 1; }
                    cnt = -1;

                }
            }
        }

        if (static_cast<int>(ncol) == 3) { ncol = column[0]; cnt = 0; }

        if (static_cast<int>(ncol) == 2)
        {
            break;
        }

        cnt++;
        num++;
    }
    //history,nameに使われているメモリを解放
    DumpFile.close();
    history.clear();
    std::map<int, EventInfo>(history).swap(history);
    name.clear();
    name.shrink_to_fit();

    std::cout<<"Finished\nConverting to Pulse\n";

    // 出力ディレクトリを作成
    std::filesystem::create_directories(DataPath+"/PulseCPP/Ch0");
    std::filesystem::create_directories(DataPath+"/PulseCPP/CH1");

	nlohmann::ordered_json json_obj;

    int n_abs = InputJson["n_abs"];
        
    const std::vector<double> Block = linspace(-1, 1, n_abs+1);

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
            for(size_t i=0;i< inner_pair.second.E_deposit.size();i++)
            {
                int Pixel = InBlock(Block, inner_pair.second.x_deposit[i], inner_pair.second.y_deposit[i], inner_pair.second.z_deposit[i]);
                BlockDeposit[Pixel - 1] += inner_pair.second.E_deposit[i];
            }

            for(auto& block_deposit:BlockDeposit)
            {
                InputPara.positions.push_back(block_deposit * 1e6) ;
            }
        }
        Eigen::MatrixXd pulse_total_0;
        Eigen::MatrixXd pulse_total_1;
        std::cout << "Before Model\n";
        Model(InputPara, pulse_total_0, pulse_total_1);
        std::cout << "AfterModel\n";
        Eigen::VectorXd pulse_0 = pulse_total_0.rowwise().sum();
        Eigen::VectorXd pulse_1 = pulse_total_0.rowwise().sum();

        std::string ChFile_0 = DataPath + "/PulseCPP/Ch0/CH0_" + std::to_string(outer_pair.first) + ".dat";
        std::ofstream outFile_0(ChFile_0);
        if (!outFile_0) {
            std::cerr << "ファイルを開けませんでした。" << std::endl;
            return;
        }
        for (int i = 0; i < pulse_0.size(); ++i) {
            outFile_0 << pulse_0(i) << std::endl;
        }
        outFile_0.close();

        std::string ChFile_1 = DataPath + "/PulseCPP/Ch1/CH1_" + std::to_string(outer_pair.first) + ".dat";
        std::ofstream outFile_1(ChFile_0);
        if (!outFile_1) {
            std::cerr << "ファイルを開けませんでした。" << std::endl;
            return;
        }
        for (int i = 0; i < pulse_0.size(); ++i) {
            outFile_1 << pulse_1(i) << std::endl;
        }
        outFile_1.close();
        Counter++;
    }

	std::cout << "Completed!\n";

    return;
}
