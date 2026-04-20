#ifndef _PARSE_TOML_H_
#define _PARSE_TOML_H_


struct TomlValue {
	std::string str;
	double number = 0;
	bool isString = false;
	bool isNumber = false;

	std::vector<TomlValue> list;
	std::map<std::string, TomlValue> table;
};


class SimpleToml {
public:
	std::map<std::string, TomlValue> root;

	bool load(const std::string& filename);


//private:
	TomlValue parseValue(const std::string& v);
	void setValue(
		const std::vector<std::string>& path,
		const std::string& key,
		const TomlValue& val);

	static void trim(std::string& s);
	static bool isNumberLiteral(const std::string& s);
	static void split(
		const std::string& s, char delim, 
		std::vector<std::string>& out);
	static void splitList(
		const std::string& s, 
		std::vector<std::string>& out);
};

#endif
