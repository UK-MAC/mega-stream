
CC = icc
CFLAGS = -std=c99 -O3

FTN = ifort
FFLAGS = -O3

OMP = -qopenmp

default: mega-stream

mega-stream: mega-stream.c
	$(CC) $(CFLAGS) $(OMP) $^ -o $@

mega-stream-omp4: mega-stream-omp4.c
	$(CC) $(CFLAGS) $(OMP) $^ -o $@

mega-stream-ftn: mega-stream.f90
	$(FTN) $(FFLAGS) $(OMP) $^ -o $@

.PHONY: clean

clean:
	rm -f mega-stream
