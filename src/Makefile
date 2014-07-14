# Simple Makefile for Kvazaar HEVC encoder

ifeq (, $(ARCH))
  ARCH = $(shell uname -m)
endif
SYSTEM = $(shell uname -s)
ASMFLAGS =
DFLAGS = -O2 -g -Werror -ftree-vectorize

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

# Do not use variable length arrays because they don't work in Visual Studio 2013.
WARNINGS = -Wall -Wtype-limits -Wvla

CC = gcc
CCFLAGS = $(DFLAGS) -std=gnu99 -I. -I./strategies $(WARNINGS)
LDFLAGS += -lm
LD = gcc -pthread -lrt
OBJS = interface_main.o \
  encmain.o \
  bitstream.o \
  cabac.o \
  checkpoint.o \
  config.o \
  context.o \
  cu.o \
  encoder.o \
  encoderstate.o \
  filter.o \
  inter.o \
  intra.o \
  nal.o \
  imagelist.o \
  rdo.o \
  sao.o \
  scalinglist.o \
  search.o \
  strategyselector.o \
  tables.o \
  threadqueue.o \
  transform.o \
  encoder_state-bitstream.o \
  encoder_state-ctors_dtors.o \
  encoder_state-geometry.o \
  image.o \
  videoframe.o \
  strategies/strategies-picture.o \
  strategies/strategies-nal.o \
  strategies/generic/nal-generic.o \
  strategies/generic/picture-generic.o \
  strategies/sse2/picture-sse2.o \
  strategies/sse41/picture-sse41.o \
  strategies/altivec/picture-altivec.o \
  strategies/avx2/picture-avx2.o
PROG  = ./kvazaar
PROGS = $(PROG)

DEPS = $(OBJS:.o=.d)


all: $(PROGS)
.PHONY: all clean


# Compile files in strategy directories with appropriate flags.
EXTRA_FLAGS =
ifeq ($(ARCH), ppc64)
  strategies/altivec/%.o: EXTRA_FLAGS += -maltivec
else
  strategies/sse2/%.o: EXTRA_FLAGS += -msse2
  strategies/sse41/%.o: EXTRA_FLAGS += -msse4.1
  strategies/avx2/%.o: EXTRA_FLAGS += -mavx2
endif


$(PROG): $(OBJS) $(ASMOBJS)
	$(LD) $^ $(LDFLAGS) -o $@

%.o: %.c Makefile
	$(CC) $(CCFLAGS) $(EXTRA_FLAGS) -c $< -o $@

%.d: %.c Makefile
	$(CC) $(CCFLAGS) -MF"$@" -MG -MM -MP -MT"$@" -MT"$(<:.c=.o)" "$<"


clean:
	echo $(ARCH)
	rm -f $(OBJS) $(PROGS) $(ASMOBJS) $(DEPS)

-include $(DEPS)
