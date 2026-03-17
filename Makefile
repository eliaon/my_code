CXX = g++

PYTHON_INC = $(shell python3-config --includes)
PYTHON_LIB = $(shell python3-config --embed --ldflags)

NUMPY_INC = $(shell python3 -c "import numpy; print(numpy.get_include())")

CXXFLAGS = -O3 -std=c++20 -fopenmp -Wno-deprecated-declarations \
           $(shell python3-config --includes) \
           -I$(NUMPY_INC)

SRCS = dipoleamplitude.cpp \
       dglap_cpp/AlphaStrong.cpp \
       dglap_cpp/EvolutionLO_nocoupling.cpp \
       main_parallel.cpp \
       integration.cpp \
       plot.cpp \
       utils.cpp \
       ctes.cpp

TARGET = main

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET) $(PYTHON_LIB)

clean:
	rm -f $(TARGET)
