cmake_minimum_required(VERSION 3.20)
project(code)

set(CMAKE_CXX_STANDARD 17)

# Caminho para o Boost (onde está a pasta boost/)
set(BOOST_INCLUDEDIR "C:/local/boost_1_89_0")

add_executable(code sigma.cpp)

# Incluir headers do Boost
target_include_directories(code PRIVATE ${BOOST_INCLUDEDIR})
