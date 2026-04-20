#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>

struct Config {

    explicit Config(const std::string& filename) {
        std::ifstream f(filename);
        if (!f.is_open())
            throw std::runtime_error("Config: cannot open '" + filename + "'");

        std::string line;
        while (std::getline(f, line)) {
            auto comment = line.find('#');
            if (comment != std::string::npos)
                line = line.substr(0, comment);

            auto eq = line.find('=');
            if (eq == std::string::npos) continue;

            std::string key   = trim(line.substr(0, eq));
            std::string value = trim(line.substr(eq + 1));

            if (!key.empty() && !value.empty())
                data_[key] = value;
        }
    }

    double get_double(const std::string& key, double def = 0.0) const {
        auto it = data_.find(key);
        if (it == data_.end()) {
            std::cerr << "Config: key '" << key << "' not found, using default " << def << "\n";
            return def;
        }
        return std::stod(it->second);
    }

    int get_int(const std::string& key, int def = 0) const {
        auto it = data_.find(key);
        if (it == data_.end()) {
            std::cerr << "Config: key '" << key << "' not found, using default " << def << "\n";
            return def;
        }
        return std::stoi(it->second);
    }

    std::string get_string(const std::string& key, const std::string& def = "") const {
        auto it = data_.find(key);
        if (it == data_.end()) {
            std::cerr << "Config: key '" << key << "' not found, using default '" << def << "'\n";
            return def;
        }
        return it->second;
    }

    void print() const {
        std::cout << "--- Config ---\n";
        for (const auto& [k, v] : data_)
            std::cout << "  " << k << " = " << v << "\n";
        std::cout << "---------------\n";
    }

private:
    std::map<std::string, std::string> data_;

    static std::string trim(const std::string& s) {
        size_t a = s.find_first_not_of(" \t\r\n");
        size_t b = s.find_last_not_of(" \t\r\n");
        return (a == std::string::npos) ? "" : s.substr(a, b - a + 1);
    }
};

#endif
