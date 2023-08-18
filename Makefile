HEADER = -I/usr/local/include
LIBB = -L/usr/local/lib
LIBRA = -lm -lfftw3
SOURCES = RV_recover.c
CFLAGS =

all:
	gcc -std=c99 $(CFLAGS) $(SOURCES) $(HEADER) $(LIBB) $(LIBRA) -o RV_recover
clean:
	rm -rf *.o RV_recover
