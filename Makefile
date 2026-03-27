CXX = g++

# Python
PYTHON_INC = $(shell python3-config --includes)
PYTHON_LIB = $(shell python3-config --embed --ldflags)

# NumPy
NUMPY_INC = $(shell python3 -c "import numpy; print(numpy.get_include())")

# ROOT
ROOT_INC  = $(shell root-config --cflags)
ROOT_LIB  = $(shell root-config --libs)

CXXFLAGS = -O3 -std=c++20 -fopenmp -Wno-deprecated-declarations \
           $(ROOT_INC) \
           $(PYTHON_INC) \
           -I$(NUMPY_INC)

LDFLAGS = $(PYTHON_LIB) $(ROOT_LIB) -lMinuit2

SRCS = dipoleamplitude.cpp \
       dglap_cpp/AlphaStrong.cpp \
       dglap_cpp/EvolutionLO_nocoupling.cpp \
       main_parallel.cpp \
       integration.cpp \
       plot.cpp \
       utils.cpp \
       ctes.cpp \
       GBW.cpp \
       bCGC.cpp \
       ipsat.cpp \
       correcs.cpp \
       wavefunctions.cpp

TARGET = main

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET)
