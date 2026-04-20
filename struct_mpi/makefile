CXX = mpic++
CXXFLAGS = -g -Wall -std=c++17 -I lib -I src
LDFLAGS =

SRCDIR = src
OBJDIR = obj
TARGET = test

# Рекурсивный поиск всех .cpp
SOURCES = $(shell find $(SRCDIR) -name "*.cpp")

# Преобразование src/...cpp → obj/...o
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

pict: $(TARGET)
	./$(TARGET)
	python3 picture.py

clean:
	rm -rf $(OBJDIR) $(TARGET)

.PHONY: all clean run pict