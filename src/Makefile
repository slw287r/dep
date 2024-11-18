ARCH := $(shell arch)
ifeq ($(ARCH),x86_64)
CFLAGS=-Wall -O3 -static
LIBS=-lhts -lpthread -ldeflate -lisal -lcairo -lfontconfig -lpixman-1 -lfreetype -lexpat -lzip -lbz2 -llzma -lpng -lm -ldl -lz
else
CFLAGS=-Wall -O3 -march=native
LIBS=-lhts -lisal -lm -lcairo
endif

CC=cc

SRC = $(wildcard *.c)
OBJS = $(patsubst %.c, %.o, $(SRC))

bamdep: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

test: bamdep
	./$^ -i ../dat/573.bam -o ../dat/573.png -s 'Klebsiella pneumoniae'

clean:
	rm -f bamdep *.o ../dat/573.png
