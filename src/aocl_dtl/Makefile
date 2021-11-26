
#
# Makefile
#
# 
#

.PHONY: generate clean

generate:
	$(CC) -fPIC -c *.h *.c
	$(CC) -shared -o libaocldtl.so *.o
	ar -rc libaocldtl.a *.o

clean: cleanx
cleanx:
	rm *.o *.gch *.a *.so
