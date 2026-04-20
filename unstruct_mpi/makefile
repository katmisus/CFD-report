MPICXX   = mpicxx
CXXFLAGS = -std=c++17 -O2
I        = -I include

MPI = cfd_mpi
MPI_OBJS = \
	main_mpi.o \
	io/vtk.o \
	mesh/mesh_io.o \
	mesh/build_faces.o \
	mpi/decomposition.o \
	mpi/mpi_mesh.o \
	physics/init.o \
	physics/bc.o \
	solver/solver.o \
	solver/hllc.o \
	solver/exact_riemann.o

$(MPI): $(MPI_OBJS)
	mkdir -p output
	$(MPICXX) $(CXXFLAGS) $^ -o $@

%.o: src/%.cpp
	$(MPICXX) $(CXXFLAGS) $(I) -c $< -o $@

all: $(MPI)

clean:
	rm -f $(MPI_OBJS) $(MPI)

clean_output:
	rm -rf output/*

.PHONY: all clean clean_output