CXX = g++

# Python
PYTHON_INC = $(shell python3-config --includes)
PYTHON_LIB = $(shell python3-config --embed --ldflags)

# NumPy
NUMPY_INC = $(shell python3 -c "import numpy; print(numpy.get_include())")

# ROOT
#ROOT_INC  = $(shell root-config --cflags) # liga se for usar ROOT
#ROOT_LIB  = $(shell root-config --libs) # liga se for usar ROOT

#LHAPDF
LHAPDF_CFLAGS = $(shell lhapdf-config --cflags)
LHAPDF_LIBS   = $(shell lhapdf-config --libs)

CXXFLAGS = -O3 -std=c++20 -fopenmp -Wno-deprecated-declarations \
           $(PYTHON_INC) \
           -I$(NUMPY_INC) \
           $(LHAPDF_CFLAGS)
#           $(ROOT_INC)

LDFLAGS = $(PYTHON_LIB) $(ROOT_LIB)  $(LHAPDF_LIBS)

SRCS = dipoleamplitude.cpp \
       LHAPDF.cpp \
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
       wavefunctions.cpp \
       DGLAP.cpp \
       ResolveDGLAP.cpp \
       IIM.cpp \
#       minuit.cpp # liga se for fazer ajuste de parâmetros

TARGET = main

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET)
