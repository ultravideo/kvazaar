# Simple Makefile for Kvazaar HEVC encoder

ifeq (, $(ARCH))
  ARCH = $(shell uname -m)
endif
SYSTEM = $(shell uname -s)
ASMFLAGS =
DFLAGS = -O2 -g -Werror -march=native

# ARCH related flags
ifeq ($(ARCH), x86_64)
    TARGET_CPU_BITS := 64
    TARGET_CPU_ARCH := x86
else ifeq ($(ARCH), ppc64)
    TARGET_CPU_BITS := 64
    TARGET_CPU_ARCH := ppc
else
    #safe (?) defaults
    TARGET_CPU_BITS := 32
    TARGET_CPU_ARCH := x86
endif

DFLAGS += -m$(TARGET_CPU_BITS)
LDFLAGS += -m$(TARGET_CPU_BITS)

# Windows (cygwin/mingw) specific flags
ifneq ( ,$(findstring Windows, $(OS)))
  ifeq ($(ARCH), x86_64)
    ASMFLAGS += -f win64
  else
    ASMFLAGS += -f win32
    ASMFLAGS += -DPREFIX
  endif
  DFLAGS += -D__USE_MINGW_ANSI_STDIO=1
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
CCFLAGS = $(DFLAGS) -I. -Wall -Wtype-limits
LDFLAGS += -lm
LD = gcc -pthread -lrt
OBJS = interface_main.o encmain.o bitstream.o cabac.o config.o context.o encoder.o encoderstate.o filter.o inter.o intra.o nal.o picture.o rdo.o sao.o scalinglist.o search.o strategyselector.o tables.o threadqueue.o transform.o
PROG  = ./kvazaar
PROGS = $(PROG)

DEPS = $(OBJS:.o=.d)

all: $(PROGS)
.PHONY: all clean

$(PROG): $(OBJS) $(ASMOBJS)
	$(LD) $^ $(LDFLAGS) -o $@

%.o: %.c Makefile
	$(CC) $(CCFLAGS) -c $< -o $@

%.d: %.c Makefile
	$(CC) $(CCFLAGS) -MF"$@" -MG -MM -MP -MT"$@" -MT"$(<:.c=.o)" "$<"


clean:
	echo $(ARCH)
	rm -f $(OBJS) $(PROGS) $(ASMOBJS) $(DEPS)

-include $(DEPS)
