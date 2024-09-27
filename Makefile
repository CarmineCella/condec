# Compiler
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2

# Source files for different targets
CONDEC_SRC = condec.cpp
FEATURES_SRC = features.cpp

# Object files
CONDEC_OBJ = $(CONDEC_SRC:.cpp=.o)
FEATURES_OBJ = $(FEATURES_SRC:.cpp=.o)

# Output executable names
CONDEC_TARGET = condec
FEATURES_TARGET = features

# Default target: build both executables
all: $(CONDEC_TARGET) $(FEATURES_TARGET)

# Compile condec executable
$(CONDEC_TARGET): $(CONDEC_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile features executable
$(FEATURES_TARGET): $(FEATURES_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile individual source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up object files and executables
clean:
	rm -f $(CONDEC_OBJ) $(FEATURES_OBJ) $(CONDEC_TARGET) $(FEATURES_TARGET)

.PHONY: all clean
