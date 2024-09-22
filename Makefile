# Makefile for compiling all .cpp files in the current folder

# Compiler
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2

# Find all .cpp files in the current directory
SRCS = $(wildcard *.cpp)

# Generate object file names by replacing .cpp with .o
OBJS = $(SRCS:.cpp=.o)

# Output executable name
TARGET = condec

# Default target
all: $(TARGET)

# Link the object files to create the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile the source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up object files and executable
clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
