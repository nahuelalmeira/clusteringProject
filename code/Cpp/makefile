SRC_PATH=./source
HEADER_PATH=./headers
OUT_PATH=../bin
CC=g++
CFLAGS=-Wall -o2
FILES=annealing.cpp basicFunctions.cpp graphFunctions.cpp ran2.cpp
SOURCES=$(FILES:%.cpp=$(SRC_PATH)/%.cpp)

all :
	mkdir -p ${OUT_PATH}
	${CC} --std=c++11 ${CFLAGS} -I${HEADER_PATH} ${SOURCES} -o annealing 
