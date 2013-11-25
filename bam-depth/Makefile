#
# Makefile
# uqcskenn, 2013-11-22 00:01
#
CC = gcc
CFLAGS = -g -O2 -Wall -pthread
HTS = $(HOME)/local/htslib
LIBS = -lz -lpthread -lm
EXECUTABLE = bamdepth

SOURCES = main.c bedidx.c

install: $(EXECUTABLE)
	ln -s bam-depth/$(EXECUTABLE) ..

$(EXECUTABLE): $(SOURCES)
	$(CC) $(CFLAGS) -L$(HTS) -I$(HTS) -o $@ $^ $(HTS)/libhts.a $(LIBS)

all: $(EXECUTABLE)


# vim:ft=make
#

