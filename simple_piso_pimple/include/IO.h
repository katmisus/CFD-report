#pragma once
#include "Types.h"
#include <string>


void saveFieldsVTK(const Fields& f,
                   const SimParams& par,
                   const std::string& dir,
                   int step,
                   double t);


void saveFieldsCSV(const Fields& f,
                   const SimParams& par,
                   const std::string& dir,
                   int step,
                   double t);
