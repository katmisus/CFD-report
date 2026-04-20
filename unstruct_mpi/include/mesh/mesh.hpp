#ifndef MESH_HPP
#define MESH_HPP

#include <vector>

// CHECK: MESH_STRUCTS

struct Vec4 {
    double v[4] = {0,0,0,0};
    double&       operator[](int i)       { return v[i]; }
    const double& operator[](int i) const { return v[i]; }
    void fill(double val) { for(int i=0;i<4;i++) v[i]=val; }
};

inline Vec4 operator+(const Vec4& a, const Vec4& b){ Vec4 r; for(int i=0;i<4;i++) r.v[i]=a.v[i]+b.v[i]; return r; }
inline Vec4 operator-(const Vec4& a, const Vec4& b){ Vec4 r; for(int i=0;i<4;i++) r.v[i]=a.v[i]-b.v[i]; return r; }
inline Vec4 operator*(const Vec4& a, double s)     { Vec4 r; for(int i=0;i<4;i++) r.v[i]=a.v[i]*s;      return r; }
inline Vec4 operator*(double s, const Vec4& a)     { return a*s; }
inline Vec4& operator+=(Vec4& a, const Vec4& b)    { for(int i=0;i<4;i++) a.v[i]+=b.v[i]; return a; }
inline Vec4& operator-=(Vec4& a, const Vec4& b)    { for(int i=0;i<4;i++) a.v[i]-=b.v[i]; return a; }

struct Node {
    double x, y;
};

struct Face {
    int left  = -1;   
    int right = -1;   
    int node0 = -1;
    int node1 = -1;

    double nx = 0, ny = 0;

    double length = 0;

    double mx = 0, my = 0;

    enum class BC {
        Interior, Wall, Inflow, Outflow, Farfield, Symmetry,
        MPIBound   // right
    } bc = BC::Interior;

    // Для MPIBound: ранг процесса-соседа
    int remote_rank = -1;

    bool is_boundary()  const { return right < 0; }
    bool is_mpi_bound() const { return bc == BC::MPIBound; }
};


struct Cell {
    std::vector<int> node_ids;  
    std::vector<int> face_ids;  

    double vol    = 0;      
    double cx = 0, cy = 0; 

    Vec4 U;    
    Vec4 res; 
};

struct Mesh {
    std::vector<Node> nodes;
    std::vector<Cell> cells;
    std::vector<Face> faces;

    int nc() const { return (int)cells.size(); }
    int nf() const { return (int)faces.size(); }
    int nn() const { return (int)nodes.size(); }

    void zero_res() { for(auto& c: cells) c.res.fill(0.0); }
};

#endif
