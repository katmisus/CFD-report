CXX      = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra -I include
TARGET   = inco_solver

SRCDIR = src
OBJDIR = obj

SOURCES  = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS  = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR) $(TARGET)

