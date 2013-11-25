#
# Makefile
# uqcskenn, 2013-11-22 00:01
#
CC = gcc
CFLAGS = -g -O2 -Wall -pthread
HTS = $(HOME)/local/htslib
LIBS = -lz -lpthread -lm
EXECUTABLE = bam-depth

all:
	$(CC) $(CFLAGS) -L$(HTS) -I$(HTS) -o $(EXECUTABLE) main.c bedidx.c $(HTS)/libhts.a $(LIBS)


# vim:ft=make
#

