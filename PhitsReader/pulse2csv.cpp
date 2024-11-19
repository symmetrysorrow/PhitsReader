#pragma once
#include"pulse2csv.h"

void add_index_range(vectori& indices, int beg, int end, int inc)
{
    for (int i = beg; i <= end; i += inc)
    {
        indices.push_back(i);
    }
}

void add_index_const(vectori& indices, int value, size_t numel)
{
    while (numel--)
    {
        indices.push_back(value);
    }
}

void append_vector(vectord& vec, const vectord& tail)
{
    vec.insert(vec.end(), tail.begin(), tail.end());
}

vectord subvector_reverse(const vectord& vec, int idx_end, int idx_start)
{
    vectord result(&vec[idx_start], &vec[idx_end + 1]);
    std::reverse(result.begin(), result.end());
    return result;
}

inline int max_val(const vectori& vec)
{
    return std::max_element(vec.begin(), vec.end())[0];
}

void filter(vectord B, vectord A, const vectord& X, vectord& Y, vectord& Zi)
{
    if (A.empty())
    {
        throw std::domain_error("The feedback filter coefficients are empty.");
    }
    if (std::all_of(A.begin(), A.end(), [](double coef) { return coef == 0; }))
    {
        throw std::domain_error("At least one of the feedback filter coefficients has to be non-zero.");
    }
    if (A[0] == 0)
    {
        throw std::domain_error("First feedback coefficient has to be non-zero.");
    }

    // Normalize feedback coefficients if a[0] != 1;
    auto a0 = A[0];
    if (a0 != 1.0)
    {
        std::transform(A.begin(), A.end(), A.begin(), [a0](double v) { return v / a0; });
        std::transform(B.begin(), B.end(), B.begin(), [a0](double v) { return v / a0; });
    }

    size_t input_size = X.size();
    size_t filter_order = std::max(A.size(), B.size());
    B.resize(filter_order, 0);
    A.resize(filter_order, 0);
    Zi.resize(filter_order, 0);
    Y.resize(input_size);

    const double* x = &X[0];
    const double* b = &B[0];
    const double* a = &A[0];
    double* z = &Zi[0];
    double* y = &Y[0];

    for (size_t i = 0; i < input_size; ++i)
    {
        size_t order = filter_order - 1;
        while (order)
        {
            if (i >= order)
            {
                z[order - 1] = b[order] * x[i - order] - a[order] * y[i - order] + z[order];
            }
            --order;
        }
        y[i] = b[0] * x[i] + z[0];
    }
    Zi.resize(filter_order - 1);
}

void filtfilt(vectord B, vectord A, const vectord& X, vectord& Y)
{
    using namespace Eigen;

    int len = X.size();     // length of input
    int na = A.size();
    int nb = B.size();
    int nfilt = (nb > na) ? nb : na;
    int nfact = 3 * (nfilt - 1); // length of edge transients

    if (len <= nfact)
    {
        throw std::domain_error("Input data too short! Data must have length more than 3 times filter order.");
    }

    // set up filter's initial conditions to remove DC offset problems at the
    // beginning and end of the sequence
    B.resize(nfilt, 0);
    A.resize(nfilt, 0);

    vectori rows, cols;
    //rows = [1:nfilt-1           2:nfilt-1             1:nfilt-2];
    add_index_range(rows, 0, nfilt - 2);
    if (nfilt > 2)
    {
        add_index_range(rows, 1, nfilt - 2);
        add_index_range(rows, 0, nfilt - 3);
    }
    //cols = [ones(1,nfilt-1)         2:nfilt-1          2:nfilt-1];
    add_index_const(cols, 0, nfilt - 1);
    if (nfilt > 2)
    {
        add_index_range(cols, 1, nfilt - 2);
        add_index_range(cols, 1, nfilt - 2);
    }
    // data = [1+a(2)         a(3:nfilt)        ones(1,nfilt-2)    -ones(1,nfilt-2)];

    auto klen = rows.size();
    vectord data;
    data.resize(klen);
    data[0] = 1 + A[1];  int j = 1;
    if (nfilt > 2)
    {
        for (int i = 2; i < nfilt; i++)
            data[j++] = A[i];
        for (int i = 0; i < nfilt - 2; i++)
            data[j++] = 1.0;
        for (int i = 0; i < nfilt - 2; i++)
            data[j++] = -1.0;
    }

    vectord leftpad = subvector_reverse(X, nfact, 1);
    double _2x0 = 2 * X[0];
    std::transform(leftpad.begin(), leftpad.end(), leftpad.begin(), [_2x0](double val) {return _2x0 - val; });

    vectord rightpad = subvector_reverse(X, len - 2, len - nfact - 1);
    double _2xl = 2 * X[len - 1];
    std::transform(rightpad.begin(), rightpad.end(), rightpad.begin(), [_2xl](double val) {return _2xl - val; });

    double y0;
    vectord signal1, signal2, zi;

    signal1.reserve(leftpad.size() + X.size() + rightpad.size());
    append_vector(signal1, leftpad);
    append_vector(signal1, X);
    append_vector(signal1, rightpad);

    // Calculate initial conditions
    MatrixXd sp = MatrixXd::Zero(max_val(rows) + 1, max_val(cols) + 1);
    for (size_t k = 0; k < klen; ++k)
    {
        sp(rows[k], cols[k]) = data[k];
    }
    auto bb = VectorXd::Map(B.data(), B.size());
    auto aa = VectorXd::Map(A.data(), A.size());
    MatrixXd zzi = (sp.inverse() * (bb.segment(1, nfilt - 1) - (bb(0) * aa.segment(1, nfilt - 1))));
    zi.resize(zzi.size());

    // Do the forward and backward filtering
    y0 = signal1[0];
    std::transform(zzi.data(), zzi.data() + zzi.size(), zi.begin(), [y0](double val) { return val * y0; });
    filter(B, A, signal1, signal2, zi);
    std::reverse(signal2.begin(), signal2.end());
    y0 = signal2[0];
    std::transform(zzi.data(), zzi.data() + zzi.size(), zi.begin(), [y0](double val) { return val * y0; });
    filter(B, A, signal2, signal1, zi);
    Y = subvector_reverse(signal1, signal1.size() - nfact - 1, nfact);
}

std::tuple<double, double, int> peak_c(const std::vector<double>& data, int presamples, int w_max, int x_av, int w_av) {
    // w_max �͈͓̔�ōő�l�Ƃ��̃C���f�b�N�X�������
    auto start_iter = data.begin() + presamples;
    auto end_iter = data.begin() + presamples + w_max;

    auto peak_iter = std::max_element(start_iter, end_iter);
    double peak = *peak_iter;  // �ő�l
    int peak_index = std::distance(data.begin(), peak_iter);  // �ő�l�̃C���f�b�N�X

    // ���ϒl��v�Z (peak_index - x_av ���� w_av �͈̔�)
    int av_start = peak_index - x_av;
    double peak_av = std::accumulate(data.begin() + av_start, data.begin() + av_start + w_av, 0.0) / w_av;

    // ���ʂ�Ԃ� (�s�[�N�l�A���ӂ̕��ϒl�A�s�[�N�̃C���f�b�N�X)
    return { peak, peak_av, peak_index };
}

Eigen::VectorXd coder_array_to_eigen(const coder::array<double, 2U>& coder_arr) {
    // coder_arr �̃T�C�Y��m�F
    if (coder_arr.size(0) == 0) {
        throw std::invalid_argument("Input coder::array is empty.");
    }

    // Eigen::VectorXd �������
    Eigen::VectorXd eigen_vec(coder_arr.size(0));

    // coder::array ���� Eigen �x�N�g���ւ̃f�[�^�R�s�[
    for (int i = 0; i < coder_arr.size(0); ++i) {
        eigen_vec(i) = coder_arr[i]; // �񎟌��z��̂��ߗ�C���f�b�N�X��0
    }

    return eigen_vec;
}

std::tuple<double, int, int> risetime(const std::vector<double>& data, double peak, int peak_index, double rise_high, double rise_low, double rate) {
    int rise_90 = 0;
    int rise_10 = 0;

    // peak_index ����t�����ɃX�L�������Arise_high �ɑ�������ʒu��T��
    for (int i = peak_index; i >= 0; --i) {
        if (data[i] <= peak * rise_high) {
            rise_90 = i;
            break;
        }
    }

    // rise_90 ����t�����ɃX�L�������Arise_low �ɑ�������ʒu��T��
    for (int j = rise_90; j >= 0; --j) {
        if (data[j] <= peak * rise_low) {
            rise_10 = j;
            break;
        }
    }

    // ���C�Y�^�C����v�Z
    double rise = (rise_90 - rise_10) / rate;

    // ���ʂ�Ԃ� (���C�Y�^�C���Arise_10�Arise_90 �̃C���f�b�N�X)
    return { rise, rise_10, rise_90 };
}

// Eigen::VectorXd �� coder::array<double, 2U> �ɕϊ�����֐�
void eigen_to_coder_array(const Eigen::VectorXd& eigen_vec, coder::array<double, 2U>& coder_arr) {
    // coder::array �̃T�C�Y��ݒ�
    coder_arr.set_size(eigen_vec.size(), 1);  // (�s��, ��)

    // Eigen �x�N�g������ coder::array �ւ̃f�[�^�R�s�[
    for (int i = 0; i < eigen_vec.size(); ++i) {
        coder_arr[i + 1] = eigen_vec(i); // �C���f�b�N�X��1�n�܂�ɒ���
    }
}

std::pair<std::vector<double>, std::vector<double>> MakeCoeff(const InputParameters& InputPara) {
    int N = 2;

    coder::array<double, 2U> bd;
    coder::array<double, 2U> ad;
    bd.set_size(1, 3);  // �s��̃T�C�Y��ݒ� (�s, ��)
    ad.set_size(1, 3);  // �s��̃T�C�Y��ݒ� (�s, ��)

    // besselfi �֐��̌Ăяo��
    besselfi(N, InputPara.cutoff, InputPara.rate, bd, ad);

    vectord b_coeff = { bd[0],bd[1],bd[2] };
    vectord a_coeff = { ad[0],ad[1],ad[2] };

    return { b_coeff,a_coeff };
}

std::vector<double> ApplyFilter(const Eigen::VectorXd& Vector,const std::vector<double>& b_coeff, const std::vector<double>& a_coeff)
{
    std::vector<double> signal(Vector.data(), Vector.data() + Vector.size());

    vectord y_filter_out;
    filtfilt(b_coeff, a_coeff, signal, y_filter_out);

    return y_filter_out;
}

std::tuple<int, double, double> GetPulseInfo(const int& Event, const std::vector<double>& Signal)
{
    auto [peak, peak_av, peak_index] = peak_c(Signal, 0, 10000, 10, 100);

    auto [rise_time, rise_10, rise_90] = risetime(Signal, peak_av, peak_index, 0.9, 0.1, 1000000);

    return { Event, peak_av,rise_time };
}



