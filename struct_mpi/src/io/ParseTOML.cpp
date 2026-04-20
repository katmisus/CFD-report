#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <cctype>
#include "ParseTOML.h"

bool SimpleToml::load(const std::string& filename) {
        std::ifstream f(filename);
        if (!f) return false;

        std::string line;
        std::vector<std::string> currentPath;

        while (std::getline(f, line)) {
            trim(line);

            // skip blank/ comment
            if (line.empty() || line[0] == '#')
                continue;

            // SECTION ----------
            if (line.front() == '[' && line.back() == ']') {
                std::string sec = line.substr(1, line.size() - 2);
                trim(sec);

                currentPath.clear();
                split(sec, '.', currentPath);
                continue;
            }

            // KEY = VALUE -------
            size_t eq = line.find('=');
            if (eq == std::string::npos)
                continue;

            std::string key = line.substr(0, eq);
            std::string val = line.substr(eq + 1);

            trim(key);
            trim(val);

            TomlValue v = parseValue(val);
            setValue(currentPath, key, v);
        }

        return true;
}

    // ----------------------- PARSE VALUE ------------------------

TomlValue SimpleToml::parseValue(const std::string& v) {
        TomlValue tv;

        // STRING: "text"
        if (v.size() > 1 && v.front() == '"' && v.back() == '"') {
            tv.str = v.substr(1, v.size() - 2);
            tv.isString = true;
            return tv;
        }

        // LIST: [ ... ]
        if (v.size() > 1 && v.front() == '[' && v.back() == ']') {
            std::string inside = v.substr(1, v.size() - 2);

            std::vector<std::string> parts;
            splitList(inside, parts);

            for (auto &p : parts) {
                trim(p);
                tv.list.push_back(parseValue(p));
            }
            return tv;
        }

        // NUMBER
        if (isNumberLiteral(v)) {
            tv.number = std::stod(v);
            tv.isNumber = true;
            return tv;
        }

        // RAW STRING (no quotes)
        tv.str = v;
        tv.isString = true;
        return tv;
}

    // ----------------------- SET VALUE --------------------------

void SimpleToml::setValue(const std::vector<std::string>& path,
                  const std::string& key,
                  const TomlValue& val){
        TomlValue* node = &root[path.empty() ? "" : path[0]];

        for (size_t i = 1; i < path.size(); i++) {
            node = &node->table[path[i]];
        }

        node->table[key] = val;
}

    // ----------------------- UTILS -------------------------------

void SimpleToml::trim(std::string& s) {
        while (!s.empty() && isspace(s.front())) s.erase(s.begin());
        while (!s.empty() && isspace(s.back())) s.pop_back();
    }

bool SimpleToml::isNumberLiteral(const std::string& s) {
	if (s.empty()) return false;

    size_t i = 0;

    // optional sign
    if (s[i] == '-' || s[i] == '+') {
        i++;
        if (i == s.size()) return false; // just "+" or "-"
    }

    bool digit = false;
    bool dot = false;

    for (; i < s.size(); i++) {
        char c = s[i];

        if (isdigit(c)) {
            digit = true;
        }
        else if (c == '.' && !dot) {
            dot = true;
        }
        else {
            return false;
        }
    }

    return digit;     

}

void SimpleToml::split(const std::string& s, char delim, std::vector<std::string>& out) {
        out.clear();
        std::stringstream ss(s);
        std::string item;
        while (std::getline(ss, item, delim)) {
            trim(item);
            out.push_back(item);
        }
}

void SimpleToml::splitList(const std::string& s, std::vector<std::string>& out) {
        out.clear();
        std::stringstream ss(s);
        std::string item;

        std::string buf;
        int bracket = 0;

        for (char c : s) {
            if (c == ',' && bracket == 0) {
                out.push_back(buf);
                buf.clear();
            } else {
                buf += c;
                if (c == '[') bracket++;
                if (c == ']') bracket--;
            }
        }
        if (!buf.empty())
            out.push_back(buf);
}


