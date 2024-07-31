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

#define Python

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

#ifdef Python
extern "C" __declspec(dllexport) void MakeOutput(const char* DatPath) {

    std::string path(DatPath);
#else
void main(){
    std::string path = "dumpall.dat";
#endif

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

    std::ifstream file(DatPath, std::ios::binary);

    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << path << std::endl;
        return;
    }

    std::cout << "Processing file...\n";
        
    std::string line;
    while (std::getline(file, line)) //ファイルの各行ごとに実行
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
                int iityp = static_cast<int>(ityp);
                if (iityp == 14 || iityp == 12 || iityp == 13)
                {
                    if ((energy >= emin_electron && iityp == 12) || (energy >= emin_photon && iityp == 14) || (energy >= emin_electron && iityp == 13))
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
    file.close();
    history.clear();
    std::map<int, EventInfo>(history).swap(history);
    name.clear();
    name.shrink_to_fit();

    std::cout<<"Finished\nWriting output.json...\n";

    try {
        nlohmann::ordered_json json_obj;

        for (const auto& outer_pair : batch) {
            for (const auto& inner_pair : outer_pair.second) {
                json_obj[std::to_string(outer_pair.first)][std::to_string(inner_pair.first)]["ityp"] = inner_pair.second.ityp;
                json_obj[std::to_string(outer_pair.first)][std::to_string(inner_pair.first)]["x"] = inner_pair.second.x;
                json_obj[std::to_string(outer_pair.first)][std::to_string(inner_pair.first)]["y"] = inner_pair.second.y;
                json_obj[std::to_string(outer_pair.first)][std::to_string(inner_pair.first)]["z"] = inner_pair.second.z;
                json_obj[std::to_string(outer_pair.first)][std::to_string(inner_pair.first)]["E"] = inner_pair.second.E;
                json_obj[std::to_string(outer_pair.first)][std::to_string(inner_pair.first)]["x_deposit"] = inner_pair.second.x_deposit;
                json_obj[std::to_string(outer_pair.first)][std::to_string(inner_pair.first)]["y_deposit"] = inner_pair.second.y_deposit;
                json_obj[std::to_string(outer_pair.first)][std::to_string(inner_pair.first)]["z_deposit"] = inner_pair.second.z_deposit;
                json_obj[std::to_string(outer_pair.first)][std::to_string(inner_pair.first)]["E_deposit"] = inner_pair.second.E_deposit;
            }
        }

        std::ofstream output_stream(output_file);
        output_stream << std::setw(4) << json_obj;
        output_stream.close();

        std::cout << "Completed!\n";
    }
    catch (const std::exception& e) {
        std::cerr << "Error writing to JSON file: " << e.what() << std::endl;
        return;
    }

    return;
}
