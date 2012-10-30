CXX := clang++
LDLIBS := -lboost_serialization

EXE := import_reads

all:$(EXE)

import_reads:import_reads.o BaseVecVec.o util.o 
	$(CXX) -o$@ $(LDFLAGS) $+ $(LDLIBS)

clean:
	rm -f *.o $(EXE)
