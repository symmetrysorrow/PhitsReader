#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <filesystem>
#include <unordered_set>
#include <nlohmann/json.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

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

// Vectorをptreeに変換するヘルパー関数
boost::property_tree::ptree to_ptree(const std::vector<double>& vec) {
    boost::property_tree::ptree pt;
    for (const auto& elem : vec) 
    {
        boost::property_tree::ptree item;
        item.put("", elem);
        pt.push_back(std::make_pair("", item));
    }
    return pt;
}

// EventInfoをptreeに変換する関数
boost::property_tree::ptree to_ptree(const EventInfo& event) {
    boost::property_tree::ptree pt;
    pt.put("ityp", event.ityp);
    pt.add_child("x", to_ptree(event.x));
    pt.add_child("y", to_ptree(event.y));
    pt.add_child("z", to_ptree(event.z));
    pt.add_child("E", to_ptree(event.E));
    pt.add_child("x_deposit", to_ptree(event.x_deposit));
    pt.add_child("y_deposit", to_ptree(event.y_deposit));
    pt.add_child("z_deposit", to_ptree(event.z_deposit));
    pt.add_child("E_deposit", to_ptree(event.E_deposit));
    return pt;
}
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

int main() {

    int counter = 0;

    //定数パラメーター
    constexpr double emin_electron = 0.1;
    constexpr double emin_photon = 0.001;

    //all plot[0], one plot[eventnumber], no plot[-1]
    constexpr double event_number = -1;
   
    std::string DatPath = "dumpall.dat";
    std::string output_file = "output_cpp.json";
#if 0
    if (std::filesystem::exists("input.json"))
    {
        boost::property_tree::ptree pt;
        boost::property_tree::read_json("input.json", pt);
        auto output = pt.get<std::string>("output");
        DatPath = "G:/tagawa/phits/trace/output_8/" + output + "/dumpall.dat";
    }
    else
    {
        std::cout << "Can't find input.json\n";
        return 1;
    }
#endif
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
        std::cerr << "Failed to open file: " << DatPath << std::endl;
        return 1;
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
        nlohmann::ordered_json JsonObj;

        for (const auto& outer_pair : batch) {
            for (const auto& inner_pair : outer_pair.second) {
                JsonObj[std::to_string(outer_pair.first)][std::to_string(inner_pair.first)]["ityp"] = inner_pair.second.ityp;
                JsonObj[std::to_string(outer_pair.first)][std::to_string(inner_pair.first)]["x"] = inner_pair.second.x;
                JsonObj[std::to_string(outer_pair.first)][std::to_string(inner_pair.first)]["y"] = inner_pair.second.y;
                JsonObj[std::to_string(outer_pair.first)][std::to_string(inner_pair.first)]["z"] = inner_pair.second.z;
                JsonObj[std::to_string(outer_pair.first)][std::to_string(inner_pair.first)]["E"] = inner_pair.second.E;
                JsonObj[std::to_string(outer_pair.first)][std::to_string(inner_pair.first)]["x_deposit"] = inner_pair.second.x_deposit;
                JsonObj[std::to_string(outer_pair.first)][std::to_string(inner_pair.first)]["y_deposit"] = inner_pair.second.y_deposit;
                JsonObj[std::to_string(outer_pair.first)][std::to_string(inner_pair.first)]["z_deposit"] = inner_pair.second.z_deposit;
                JsonObj[std::to_string(outer_pair.first)][std::to_string(inner_pair.first)]["E_deposit"] = inner_pair.second.E_deposit;
            }
        }

        std::ofstream OutputStream(output_file);
        OutputStream << std::setw(4) << JsonObj<< std::endl;
        OutputStream.close();

        std::cout << "Completed!\n";

        int some;
        std::cin >> some;

    }
    catch (const std::exception& e) {
        std::cerr << "Error writing to JSON file: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
