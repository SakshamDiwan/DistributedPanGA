CC = cc
CFLAGS = -O3 -Wall -I./lib
LIBS = -lpthread -lm -lz

all: distributedpanga

distributedpanga: src/main.c
	$(CC) $(CFLAGS) -o distributedpanga src/main.c $(LIBS)

clean:
	rm -f distributedpanga

