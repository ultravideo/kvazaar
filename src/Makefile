# Simple Makefile for Kvazaar HEVC encoder


BITS = $(shell getconf LONG_BIT)

# Detect 32/64bit system
ifeq ($(BITS), 64)
    ASMOBJS = test64.o
else
    ASMOBJS = test.o
endif


DFLAGS =
ifeq ($(SYSTEM),Darwin)
DFLAGS += -D__APPLE__
endif
ifeq ($(SYSTEM),FreeBSD)
DFLAGS += -D__FREEBSD__
endif
ifeq ($(SYSTEM),SunOS)
DFLAGS += -D__SOLARIS__
endif

CC = gcc
CCFLAGS = $(DFLAGS) -I. -I/usr/local/include -L/usr/local/lib
LDFLAGS = -lm
LD = gcc
YASM = yasm
OBJS = interface_main.o encmain.o bitstream.o cabac.o config.o context.o debug.o encoder.o filter.o inter.o intra.o nal.o picture.o rdo.o sao.o search.o transform.o
PROG  = ./kvazaar
PROGS = $(PROG)

$(PROG): $(OBJS) $(ASMOBJS)
	$(LD) $^ $(LDFLAGS) -o $@

test.o: x86/test.asm
	$(YASM) -f elf x86/test.asm -o test.o

test64.o: x64/test64.asm
	$(YASM) -f elf64 x64/test64.asm -o test64.o

%.o: %.c
	$(CC) $(CCFLAGS) -c $< -o $@

all: $(PROGS)
