CXX := g++
CXXFLAGS := -std=c++20 -Wall -Wextra -O2
ROOTFLAGS := $(shell root-config --cflags --libs)

SRCS := simple.C test.C
TARGET := simple

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) -o $@ $(SRCS) $(CXXFLAGS) $(ROOTFLAGS)

clean:
	rm -f $(TARGET)
