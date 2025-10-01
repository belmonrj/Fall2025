# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++20 -Wall -Wextra -O2
ROOTFLAGS := $(shell root-config --cflags --libs)

# Sources and executable
SRCS := simple.C test.C flow_functions.C
OBJS := $(SRCS:.C=.o)
TARGET := simple

# Default target
all: $(TARGET)

# Link the program
$(TARGET): $(OBJS)
	$(CXX) -o $@ $^ $(ROOTFLAGS)

# Compile object files
%.o: %.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

# Clean up
clean:
	rm -f $(TARGET) $(OBJS)
