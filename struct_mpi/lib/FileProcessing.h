#ifndef _FILE_PROCESSING_H_
#define _FILE_PROCESSING_H_

#include "Types.h"

void SaveFieldToCSV(const Field& W,
                    const std::vector<double>& x,
                    const std::vector<double>& y,
                    const double& time,
                    const std::string& filename,
                    bool append);

void SaveFluxToCSV(const Field& Flux,
                   const std::vector<double>& x,
                   const std::vector<double>& y,
                   const double& time,
                   const std::string& filename,
                   int dir);

fs::path CreateDirFromPath(const std::string& file_path);
/*void MoveToChache(const std::string& source_root_path, const std::string& cache_root_path);

void CleanDir(std::vector<std::string> folders);

void WriteToCSV(std::vector<std::vector<double>> W, std::vector<double> xc, double t, std::ofstream& file);

void ClearDirectoryContents(const std::string& target_dir_path);

void CreateDir(const std::string& root_folder_path, const std::string& new_folder_name);

void SaveAnalysisData(
    std::ofstream& file,
    double t,
    const std::vector<std::vector<double>>& W_num,
    const std::vector<double>& xc,                 
    // Параметры для точного решения:
    const std::vector<double>& W_L,
    const std::vector<double>& W_R,
    const std::vector<double>& W_star,
    double x0,
    std::vector<double>(*TrueSolveF)(std::vector<double>, std::vector<double>, std::vector<double>, double, double),
    double(*AnalisF)(std::vector<double>)
);*/


#endif
