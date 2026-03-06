CXX = g++
CXXFLAGS = -O3 -std=c++20 -fopenmp


SRCS = dipoleamplitude.cpp \
       dglap_cpp/AlphaStrong.cpp \
       dglap_cpp/EvolutionLO_nocoupling.cpp \
       main_parallel.cpp \
       integration.cpp

TARGET = main

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET)

clean:
	rm -f $(TARGET)
