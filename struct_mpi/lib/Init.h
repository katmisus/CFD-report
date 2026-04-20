#ifndef _INIT_H_
#define _INIT_H_

#include "Types.h"

void readConfig(const std::string& config_path);

void Grid(std::vector<double>& x, std::vector<double>& y,
          int offset_x, int offset_y);

void InitValues(Field& W, 
				const std::vector<double>& x, 
				const std::vector<double>& y,
				const std::string& config_path);

#endif
