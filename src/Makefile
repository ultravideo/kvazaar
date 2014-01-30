# Simple Makefile for Kvazaar HEVC encoder

ARCH = $(shell uname -p)
SYSTEM = $(shell uname -s)
ASMFLAGS =
DFLAGS =

# ARCH related flags
ifeq ($(ARCH), x86_64)
    ASMFLAGS += -DARCH_X86_64=1
else
    ASMFLAGS += -DARCH_X86_64=0
endif


# Linux specifix flags
ifeq ($(SYSTEM),Linux)
  ifeq ($(ARCH), x86_64)    
      ASMFLAGS += -f elf64
  else
      ASMFLAGS += -f elf32
  endif
endif

# OS X specific flags
ifeq ($(SYSTEM),Darwin)
  ifeq ($(ARCH), x86_64)
  ASMFLAGS += -f macho64
  else
  ASMFLAGS += -f macho32
  endif
ASMFLAGS += -DPREFIX
endif

CC = gcc
CCFLAGS = $(DFLAGS) -I. -I/usr/local/include -L/usr/local/lib
LDFLAGS = -lm
LD = gcc
YASM = yasm
ASMOBJS = cpu.o
OBJS = interface_main.o encmain.o bitstream.o cabac.o config.o context.o debug.o encoder.o filter.o inter.o intra.o nal.o picture.o rdo.o sao.o search.o transform.o
PROG  = ./kvazaar
PROGS = $(PROG)


$(PROG): $(OBJS) $(ASMOBJS)
	$(LD) $^ $(LDFLAGS) -o $@

cpu.o: x86/cpu.asm
	$(YASM) $(ASMFLAGS) x86/cpu.asm -o cpu.o

%.o: %.c
	$(CC) $(CCFLAGS) -c $< -o $@

all: $(PROGS)
