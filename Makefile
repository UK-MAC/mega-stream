
CC = icc
CFLAGS = -O3 -qopenmp

mega-stream: mega-stream.c
	$(CC) $(CFLAGS) $^ -o $@
