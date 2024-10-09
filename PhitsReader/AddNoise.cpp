#include "AddNoise.h"
#include <unsupported/Eigen/FFT>
#include <random>
#include<fstream>
#include<iostream>

Eigen::VectorXd inverseFFTNoise(const Eigen::VectorXcd& noise_spe, int noise_samples, int samples)
{
    // noise_spe のサイズが noise_samples の2倍であることを確認
    if (noise_spe.size() != 2 * noise_samples) {
        throw std::invalid_argument("noise_spe.size() must be 2 * noise_samples");
    }

    // Eigen の FFT ライブラリを使って逆フーリエ変換を実行
    Eigen::FFT<double> fft;

    // 複素数のベクトルを作成
    Eigen::VectorXcd noise_time_domain(2 * noise_samples);

    // 逆フーリエ変換を実行（周波数領域 → 時間領域）
    fft.inv(noise_time_domain, noise_spe);

    // 実部を取り出し、samples 分だけ切り詰める
    Eigen::VectorXd result(samples);
    for (int i = 0; i < samples; ++i) {
        if (i < noise_samples) {
            result[i] = 2.0 * noise_time_domain[i].real();  // 実部を取り出し2倍
        }
        else {
            result[i] = 0.0;  // 必要に応じて残りを0で埋める
        }
    }

    return result;
}

Eigen::VectorXcd random_noise(const Eigen::VectorXd& spe, unsigned seed)
{
    // 入力ベクトルを反転
    Eigen::VectorXd spe_re = spe.reverse();

    // 入力ベクトルとその反転を結合してミラーリング
    Eigen::VectorXd spe_mirror(spe.size() * 2);
    spe_mirror << spe, spe_re;

    // ランダム位相生成のための設定
    std::mt19937 generator(seed); // シードを設定した乱数生成器
    std::uniform_real_distribution<> dist(0.0, 2 * 3.141592); // 0から2πの一様分布

    // ランダム位相を使用して複素ベクトルを生成
    std::vector<std::complex<double>> complex_vec(spe_mirror.size());
    for (int i = 0; i < spe_mirror.size(); ++i) {
        double phase = dist(generator);  // ランダム位相
        complex_vec[i] = std::polar(spe_mirror[i], phase); // 複素数を作成
    }

    // ベクトルの後半部分を共役にする
    std::vector<std::complex<double>> complex_conj(spe.size());
    for (int i = 0; i < spe.size(); ++i) {
        complex_conj[i] = std::conj(complex_vec[spe.size() + i]); // 共役複素数
    }

    // 元の部分と共役部分を結合して結果を作成
    Eigen::VectorXcd result(spe.size() * 2);
    for (int i = 0; i < spe.size(); ++i) {
        result[i] = complex_vec[i]; // オリジナルの複素数部分
        result[spe.size() + i] = complex_conj[i]; // 共役部分
    }

    return result;
}

void AddNoise(const Eigen::VectorXd& Noise_dense, Eigen::VectorXd& Pulse)
{
    int NoiseSamples = Noise_dense.size();
    std::random_device rd;  // 非決定的な乱数生成
    std::mt19937 gen(rd()); // メルセンヌ・ツイスタ法による擬似乱数生成器

    // 1から10000までの範囲でランダムな整数を生成する分布を設定
    std::uniform_int_distribution<> distrib(1, 10000);
    int count= distrib(gen);

    Eigen::VectorXcd noise_spe = random_noise(Noise_dense, count) * NoiseSamples;

    Eigen::VectorXd noise = inverseFFTNoise(noise_spe, NoiseSamples, Pulse.size());
    Pulse += noise;
}

Eigen::VectorXd readLinesToEigen(const std::string& filePath)
{
    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filePath << std::endl;
        return Eigen::VectorXd();
    }

    std::vector<double> values;
    std::string line;
    while (std::getline(file, line)) {
        // 空行をスキップ
        if (line.empty()) {
            continue;
        }

        try {
            // スペースをトリムしてからdoubleに変換
            double value = std::stod(line);
            values.push_back(value);
        }
        catch (const std::invalid_argument& e) {
            std::cout << "Invalid number format: " << line << std::endl;
        }
        catch (const std::out_of_range& e) {
            std::cout << "Number out of range: " << line << std::endl;
        }
    }

    file.close();

    // VectorXdに変換
    Eigen::VectorXd result(values.size());
    for (size_t i = 0; i < values.size(); ++i) {
        result(i) = values[i];
    }

    return result;
}
