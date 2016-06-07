
CC = icc
CFLAGS = -std=c99 -O3 -qopenmp

mega-stream: mega-stream.c
	$(CC) $(CFLAGS) $^ -o $@

.PHONY: clean

clean:
	rm -f mega-stream
