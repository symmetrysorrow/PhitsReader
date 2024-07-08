#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <regex>
#include <map>
#include <filesystem>
#include <unordered_set>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

std::map<int, std::string> itype = {
    {12, "electron"},
    {13, "positron"},
    {14, "photon"}
};

inline std::string GetItype(const float& ityp) { 
    auto it = itype.find(static_cast<int>(ityp));
    if (it != itype.end()) {
        return it->second;
    }
    else {
        return "unknown";  // デフォルト値を返す
    }
};

struct EventInfo {
    int ityp;
    float x;
    float y;
    float z;
    float E;
    float x_deposit;
    float y_deposit;
    float z_deposit;
    float E_deposit;
};

boost::property_tree::ptree to_ptree(const EventInfo& e) {
    boost::property_tree::ptree pt;
    pt.put("itype", e.ityp);
    pt.put("x", e.x);
    pt.put("y", e.y);
    pt.put("z", e.z);
    pt.put("E", e.E);
    pt.put("x_deposit", e.x_deposit);
    pt.put("y_deposit", e.y_deposit);
    pt.put("z_deposit", e.z_deposit);
    pt.put("E_deposit", e.E_deposit);
    return pt;
}

std::vector<float> split_line(const std::string& line) {
    std::vector<float> column;
    std::istringstream stream(line);
    std::string token;

    while (stream >> token) {  // 空白をスキップしてトークンを取得
        try {
            column.push_back(std::stof(token));  // トークンをfloatに変換
        }
        catch (const std::invalid_argument& e) {
            std::cerr << "Invalid argument: " << token << " cannot be converted to float." << std::endl;
        }
        catch (const std::out_of_range& e) {
            std::cerr << "Out of range: " << token << " is out of range for float." << std::endl;
        }
    }

    return column;
}

int main() {

    int counter = 0;

    //parameter
    double emin_electron = 0.1;
    double emin_photon = 0.001;

    //all plot[0], one plot[eventnumber], no plot[-1]
    float event_number = -1;
    //read input.json and create values for it
    boost::property_tree::ptree pt;
    boost::property_tree::read_json("input.json", pt);
    const auto output = pt.get<std::string>("output");
    std::string path = "dumpall.dat";
    const int n_abs = pt.get<int>("n_abs");
    std::vector<double> block(n_abs + 1);

    std::map<int, EventInfo> history;
    std::map<int, std::map<int, EventInfo>> batch;
    float ncol = 1;
    std::vector<float> xyz = { 0,0,0 };
    std::vector<float> cxyz = { 0,0,0 };

    float cnt = 0;
    float num = 0;
    float nocas = 0;
    float no = 0;
    std::vector<float> reg(2);
    std::vector<float> name;
    float benergy, cenergy, ityp, nclsts, jcoll, energy_new, ncl, energy, energy_dps;
    //float iclusts;
    constexpr double start = -1.0;
    constexpr double end = 1.0;
    const double step = (end - start) / static_cast<double>(n_abs);

    for (int i = 0; i <= n_abs; ++i) {
        block[i] = start + i * step;
    }

    std::ifstream file(path, std::ios::binary);

    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << path << std::endl;
        return 1;
    }
    std::cout << "loading file...\n";
    std::vector<std::string> lines;

    std::stringstream buffer;
    buffer << file.rdbuf(); // ファイルの内容を一気に読み込む

    file.close(); // ファイルを閉じる
    std::cout << "finish loading file\nconverting to vector...\n";
    // bufferから改行ごとに行を取り出してlinesに追加
    std::string line;
    while (std::getline(buffer, line)) {
        lines.push_back(line); // 改行ごとにvectorに挿入
    }

    // bufferを解放
    buffer.str(std::string()); // stringstreamの内容をクリア
    buffer.clear();
    std::cout << "Finish converting and clear buffer\nProcessing data\n";

    for (const auto& line:lines)
    {
        counter++;
        
        if (counter % 100000 == 0)
        {
            double progress = counter == lines.size() ? 100 : static_cast<double>(counter) / lines.size() * 100.0;
            std::cout << "\rProgress: " << std::fixed << std::setprecision(0) << progress << "%";
            
            std::cout.flush();
        }
#if 0
        if (counter == 10000)
        {
            break;
        }
#endif
        std::vector<float> column=split_line(line);
#if 0
        for (const auto& col : column)
        {
            std::cout << col << "/";
        }
#endif
        //std::cout << "tag ncol " << ncol << "\n"<<"    tag cnt "<<cnt<<"\n";
        
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
                    //print_map(history);
                    //plot
                    if (no > 1)
                    {
                        batch[static_cast<int>(nocas)] = history;
                        if (event_number == nocas || static_cast<int>(event_number) == 0) {
                            //plot
                        }
                    }
                    else {
                        //plot
                    }

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
                    if (ncol == 4) { /*std::cout << "---Event" << nocas << "---\n";*/ }
                    //std::cout << num << "\n NCOL:" << ncol << "\nno:" << no << "\n ityp:" << GetItype(ityp) << "\n";
                  //  std::cout << "reg:";
                    //for (const auto& l : reg) { std::cout << l << ","; }
                    //std::cout << "\n name:";
                   // for (const auto& l : name) { std::cout << l << ","; }
                    //std::cout << "\n energy:" << benergy << "MeV \n  c-energy:" << cenergy << "MeV \n c-position:";
                  //  for (const auto& l : cxyz) { std::cout << l << ","; }
                   // std::cout << "\n -----------\n";

                    if (history.find(static_cast<int>(no)) == history.end())
                    {
                        history[static_cast<int>(no)] = { static_cast<int>(ityp),0,0,0,0,0,0,0,0 };
                    }
                    history[static_cast<int>(no)].x = cxyz[0];
                    history[static_cast<int>(no)].y = cxyz[1];
                    history[static_cast<int>(no)].z = cxyz[2];
                    history[static_cast<int>(no)].E = cenergy;

                    if (static_cast<int>(ncol) == 11)
                    {
                        history[static_cast<int>(no)].E_deposit = energy;
                        history[static_cast<int>(no)].x_deposit = cxyz[0];
                        history[static_cast<int>(no)].y_deposit = cxyz[1];
                        history[static_cast<int>(no)].z_deposit = cxyz[2];
                        //std::cout << "deposit energy:" << energy << "\n";
                       // std::cout << "--------\n";
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

                //std::cout << "nclsts:" << nclsts << "\n jcoll:" << jcoll << "\n------ \n";

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

                   // std::cout << "ityp:" << GetItype(ityp) << "\n";
                   // std::cout << "energy: " << energy << "MeV\n------";
                    if (static_cast<int>(ncl) == static_cast<int>(nclsts) - 1)
                    {
                        ncol = 13;
                        history[static_cast<int>(no)].E_deposit = energy_dps;
                        history[static_cast<int>(no)].x_deposit = cxyz[0];
                        history[static_cast<int>(no)].y_deposit = cxyz[1];
                        history[static_cast<int>(no)].z_deposit = cxyz[2];
                        //std::cout << "deposit energy:" << energy_dps << "\n------";
                    }
                    else { ncl += 1; }
                    cnt = -1;

                }
            }
        }

        if (static_cast<int>(ncol) == 3) { ncol = column[0]; cnt = 0; }

        if (static_cast<int>(ncol) == 2)
        {
            boost::property_tree::ptree pt_batch;
            for (const auto& outer_pair : batch) {
                boost::property_tree::ptree pt_inner;
                for (const auto& inner_pair : outer_pair.second) {
                    pt_inner.add_child(std::to_string(inner_pair.first), to_ptree(inner_pair.second));
                }
                pt_batch.add_child(std::to_string(outer_pair.first), pt_inner);
            }

            std::string output_file = "OutputCpp/output.json";

            std::filesystem::create_directory("OutputCpp");

            try {
                // ファイルに上書き
                boost::property_tree::write_json(output_file, pt_batch);
            }
            catch (const std::exception& e) {
                std::cerr << "Error writing to JSON file: " << e.what() << std::endl;
            }
        }

        cnt++;
        num++;
    }

    std::cout << "\nFinished!\n";
    int some;
    std::cin >> some;

    return 0;
}
