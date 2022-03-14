include make.inc
include ./GT_LPKE.inc

.SUFFIXES: .cc .c .o
.c.o:
	$(CXX)  $(CXXFLAGS) ./LIN  -I../../../../include/x86_64-unknown-linux-gnu  -c -o $@ $<
#	$(CC) $(CFLAGS) -I. -I./LIN_EIG -I../include -c -o $@ $<

.PHONY: all
all: xlapacke_test_main 

LIBRARIES =  $(LAPACKELIB) $(LAPACKLIB) $(BLASLIB) $(GTESTLIB)

xlapacke_test_main: lapacke_gtest_main.o lapacke_gtest_helper.o $(GTEST_SRC_FILE_LIST)  -ldl  $(LIBRARIES)
	$(CXX) $(CXXFLAGS)  $(LDFLAGS)  -o $@ $^ -lgfortran
#	./$@

.PHONY: clean cleanobj cleanexe
clean: cleanobj cleanexe
cleanobj:
	rm -f *.o ./LIN_EIG/*.o  ./LIN/*.o ./EIG/*.o EIG/ORTHO/*.o EIG/NEIG/*.o ./EIG/SEIG/*.o ./MISC/*.o ./LLS/*.o 
	rm -f ./LIN/DRV/*.o   ./LIN/RSEE/*.o  ./LIN/EQLB/*.o  ./SVD/*.o ./EIG/COS/*.o
cleanexe:
	rm -f x*
