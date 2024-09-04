#pragma once
#include <iostream>
#include <string>
#include <thread>
#include <atomic>
#include <chrono>
#include <vector>

class SpinProgress {
public:
    // プリセットのスピンキャラクターセット
    enum SpinType {
        Classic,
        Braille,
        Dots
    };

    // スピンキャラクターセットを選択するコンストラクタ
    SpinProgress(SpinType type = Classic)
        : is_completed(false), custom_message(""), spinner_thread(&SpinProgress::spin, this) {
        set_spin_chars(type);
    }

    ~SpinProgress() {
        if (spinner_thread.joinable()) {
            complete();
            spinner_thread.join();
        }
    }

    // スピンキャラクターセットの変更
    void set_spin_chars(SpinType type) {
        switch (type) {
        case Classic:
            spin_chars = { '|', '/', '-', '\\' };
            break;
        case Braille: {
            std::string braille_str = u8"⠈⠐⠠⢀⡀⠄⠂⠁";
            spin_chars.assign(braille_str.begin(), braille_str.end());
            break;
        }
        case Dots:
            spin_chars = { '.', 'o', 'O', '@', '*' };
            break;
        }
    }

    // 任意のメッセージを表示
    void set_message(const std::string& message) {
        custom_message = message;
    }

    // 完了を示すチェックマークとメッセージの表示
    void complete(const std::string& message = "Completed") {
        is_completed.store(true);
        if (spinner_thread.joinable()) {
            spinner_thread.join();
        }
        std::cout << "\r\033[K";  // カーソルを行の先頭に戻し、行をクリア
        std::cout << u8"✔ " << message << std::endl;
    }

private:
    std::vector<char> spin_chars;
    std::atomic<bool> is_completed;
    std::string custom_message;
    std::thread spinner_thread;

    // スピンを回転させる関数
    void spin() {
        int spin_index = 0;
        while (!is_completed.load()) {
            std::cout << "\r" << spin_chars[spin_index % spin_chars.size()] << " " << custom_message << std::flush;
            spin_index++;
            std::this_thread::sleep_for(std::chrono::milliseconds(100)); // スピンの回転速度
        }
    }
};

