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
	$^ -i ../dat/toy.bam -o ../dat/toy.bed.gz -p ../dat/toy.png -s 'Talaromyces marneffei'
	open ../dat/toy.png

clean:
	rm -f bamdep *.o ../dat/toy.bed* ../dat/toy.png ../dat/toy.svg
