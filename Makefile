
CC = icc
CFLAGS = -std=c99 -O3 -qopenmp

default: mega-stream

mega-stream: mega-stream.c
	$(CC) $(CFLAGS) $^ -o $@

mega-stream-omp4: mega-stream-omp4.c
	$(CC) $(CFLAGS) $^ -o $@

.PHONY: clean

clean:
	rm -f mega-stream
