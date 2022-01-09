
#
# Makefile
#
# 
#

.PHONY: generate clean

generate:
	gcc -fPIC -c *.h *.c
	gcc -shared -o libaocldtl.so *.o
	ar -rc libaocldtl.a *.o

clean: cleanx
cleanx:
	rm *.o *.gch *.a *.so
