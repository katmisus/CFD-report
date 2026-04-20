void DecomposeProcesses(int p,
                        int dim,
                        const std::vector<int>& global_sizes,
                        std::vector<int>& proc_dims);

int LocalSize(int N, int p_coord, int p_dim);
int Offset(int N_cells, int coord, int p_dim);