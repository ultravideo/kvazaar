# Simple Makefile for Kvazaar HEVC encoder

ifeq (, $(ARCH))
  ARCH = $(shell uname -m)
endif
SYSTEM = $(shell uname -s)
ASMFLAGS =
DFLAGS = -O2 -g -Werror

# ARCH related flags
ifeq ($(ARCH), x86_64)
    ASMFLAGS += -DARCH_X86_64=1
    DFLAGS += -m64
    LDFLAGS += -m64
else
    ASMFLAGS += -DARCH_X86_64=0
    DFLAGS += -m32
    LDFLAGS += -m32
endif

# Windows (cygwin/mingw) specific flags
ifneq ( ,$(findstring Windows, $(OS)))
  ifeq ($(ARCH), x86_64)
    ASMFLAGS += -f win64
  else
    ASMFLAGS += -f win32
    ASMFLAGS += -DPREFIX
  endif
  LDFLAGS += -lws2_32
# OS X specific flags
else ifeq ($(SYSTEM),Darwin)
  ifeq ($(ARCH), x86_64)
    ASMFLAGS += -f macho64
  else
    ASMFLAGS += -f macho32
  endif
ASMFLAGS += -DPREFIX
# Default to Linux/elf specific flags
else
  ifeq ($(ARCH), x86_64)
      ASMFLAGS += -f elf64
  else
      ASMFLAGS += -f elf32
  endif
endif

CC = gcc
CCFLAGS = $(DFLAGS) -I. -I/usr/local/include -L/usr/local/lib -Wall
LDFLAGS += -lm
LD = gcc
YASM = yasm
ASMOBJS = cpu.o
OBJS = interface_main.o encmain.o bitstream.o cabac.o config.o context.o encoder.o filter.o inter.o intra.o nal.o picture.o rdo.o sao.o scalinglist.o search.o tables.o transform.o extras/getopt.o
PROG  = ./kvazaar
PROGS = $(PROG)

DEPS = $(OBJS:.o=.d)

all: $(PROGS)
.PHONY: all clean

$(PROG): $(OBJS) $(ASMOBJS)
	$(LD) $^ $(LDFLAGS) -o $@

cpu.o: x86/cpu.asm
	$(YASM) $(ASMFLAGS) x86/cpu.asm -o cpu.o

%.o: %.c Makefile
	$(CC) $(CCFLAGS) -c $< -o $@

%.d: %.c Makefile
	$(CC) $(CCFLAGS) -MF"$@" -MG -MM -MP -MT"$@" -MT"$(<:.c=.o)" "$<"


clean:
	rm -f $(OBJS) $(PROGS) $(ASMOBJS) $(DEPS)

-include $(DEPS)
