#include "AddNoise.h"
#include <unsupported/Eigen/FFT>
#include <random>
#include<fstream>
#include<iostream>
#include<cmath>

Eigen::VectorXd inverseFFTNoise(const Eigen::VectorXcd& noise_spe, int noise_samples, int samples)
{
    // noise_spe �̃T�C�Y�� noise_samples ��2�{�ł��邱�Ƃ�m�F
    if (noise_spe.size() != 2 * noise_samples) {
        throw std::invalid_argument("noise_spe.size() must be 2 * noise_samples");
    }

    // Eigen �� FFT ���C�u������g���ċt�t�[���G�ϊ�����s
    Eigen::FFT<double> fft;

    // ���f���̃x�N�g����쐬
    Eigen::VectorXcd noise_time_domain(2 * noise_samples);

    // �t�t�[���G�ϊ�����s�i���g���̈� �� ���ԗ̈�j
    fft.inv(noise_time_domain, noise_spe);

    // ��������o���Asamples �������؂�l�߂�
    Eigen::VectorXd result(samples);
    for (int i = 0; i < samples; ++i) {
        if (i < noise_samples) {
            result[i] =  noise_time_domain[i].real();  // ��������o��2�{
        }
        else {
            result[i] = 0.0;  // �K�v�ɉ����Ďc���0�Ŗ��߂�
        }
    }

    return result;
}

Eigen::VectorXcd random_noise(const Eigen::VectorXd& spe, unsigned seed)
{
    // ���̓x�N�g���𔽓]
    Eigen::VectorXd spe_re = spe.reverse();

    // ���̓x�N�g���Ƃ��̔��]��������ă~���[�����O
    Eigen::VectorXd spe_mirror(spe.size() * 2);
    spe_mirror << spe, spe_re;

    // �����_���ʑ������̂��߂̐ݒ�
    std::mt19937 generator(seed); // �V�[�h��ݒ肵������������
    std::uniform_real_distribution<> dist(0.0, 2 * 3.141592); // 0����2�΂̈�l���z

    // �����_���ʑ���g�p���ĕ��f�x�N�g���𐶐�
    std::vector<std::complex<double>> complex_vec(spe_mirror.size());
    for (int i = 0; i < spe_mirror.size(); ++i) {
        double phase = dist(generator);  // �����_���ʑ�
        complex_vec[i] = std::polar(spe_mirror[i], phase); // ���f����쐬
    }

    // �x�N�g���̌㔼���������ɂ���
    std::vector<std::complex<double>> complex_conj(spe.size());
    for (int i = 0; i < spe.size(); ++i) {
        complex_conj[i] = std::conj(complex_vec[spe.size() + i]); // ���𕡑f��
    }

    // ���̕����Ƌ��𕔕���������Č��ʂ�쐬
    Eigen::VectorXcd result(spe.size() * 2);
    for (int i = 0; i < spe.size(); ++i) {
        result[i] = complex_vec[i]; // �I���W�i���̕��f������
        result[spe.size() + i] = complex_conj[i]; // ���𕔕�
    }

    return result;
}

void AddNoise(const Eigen::VectorXd& Noise_dense, Eigen::VectorXd& Pulse, const double& rate)
{
    int NoiseSamples = Noise_dense.size();

    double delta_f=rate/NoiseSamples;
    
    std::random_device rd;  // �񌈒�I�ȗ�������
    std::mt19937 gen(rd()); // �����Z���k�E�c�C�X�^�@�ɂ��[������������

    // 1����10000�܂ł͈̔͂Ń����_���Ȑ����𐶐����镪�z��ݒ�
    std::uniform_int_distribution<> distrib(1, 10000);
    int count= distrib(gen);

    Eigen::VectorXcd noise_spe = random_noise(Noise_dense, count);
    Eigen::VectorXcd ifft_input = noise_spe * std::sqrt(delta_f)*NoiseSamples*std::sqrt(2);

    Eigen::VectorXd noise = inverseFFTNoise(ifft_input, NoiseSamples, Pulse.size());
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
        // ��s��X�L�b�v
        if (line.empty()) {
            continue;
        }

        try {
            // �X�y�[�X��g�������Ă���double�ɕϊ�
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

    // VectorXd�ɕϊ�
    Eigen::VectorXd result(values.size());
    for (size_t i = 0; i < values.size(); ++i) {
        result(i) = values[i];
    }

    return result;
}
