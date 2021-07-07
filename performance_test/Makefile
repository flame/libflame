TOPSRCDIR = ../..
include ./make.inc
include ./GBench_Flame.inc

.SUFFIXES: .cc .c .o
.c.o:
	$(CXX)  $(CXXFLAGS) -c -o $@ $<

.cc.o:
	$(CXX)  $(CXXFLAGS) -c -o $@ $<

.PHONY: all
all: xFlame_gbench 

LIBRARIES = $(LAPACKLIB) $(BLASLIB)  $(GTESTLIB) $(GBENCHLIB)

#Linker command for AOCL linking
xFlame_gbench: Flame_gbench_main.o  $(GBENCH_SRC_FILE_LIST)  -ldl  $(LIBRARIES)
	$(CXX) $(CXXFLAGS)  $(LDFLAGS)  -o $@ $^ -lgfortran  -lpthread
#	./$@

.PHONY: clean cleanobj cleanexe
clean: cleanobj cleanexe
cleanobj:
	rm -f *.o  ./LIN/*.o
#	rm -f *.o ./LIN_EIG/*.o  ./LIN/*.o ./EIG/*.o EIG/ORTHO/*.o EIG/NEIG/*.o ./EIG/SEIG/*.o ./AUX/*.o ./LLS/*.o 
#	rm -f ./LIN/DRV/*.o   ./LIN/RSEE/*.o  ./LIN/EQLB/*.o  ./SVD/*.o ./EIG/COS/*.o
cleanexe:
	rm -f x*
