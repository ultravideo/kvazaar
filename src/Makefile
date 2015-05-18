# Simple Makefile for Kvazaar HEVC encoder

#Disable default rules
.SUFFIXES:

ifeq (, $(ARCH))
  ARCH = $(shell uname -m)
endif
SYSTEM = $(shell uname -s)
ASMFLAGS =
CFLAGS += -O2 -g -Werror -ftree-vectorize

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

CFLAGS += -m$(TARGET_CPU_BITS)
LDFLAGS += -m$(TARGET_CPU_BITS)

# Windows (cygwin/mingw) specific flags
ifneq ( ,$(findstring Windows, $(OS)))
  ifeq ($(ARCH), x86_64)
    ASMFLAGS += -f win64
    ASMFLAGS += -DHAVE_ALIGNED_STACK=1
  else
    ASMFLAGS += -f win32
    ASMFLAGS += -DPREFIX
    ASMFLAGS += -DHAVE_ALIGNED_STACK=0
  endif
  CFLAGS += -D__USE_MINGW_ANSI_STDIO=1
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
  LDFLAGS += -lrt
  ifeq ($(ARCH), x86_64)
    ASMFLAGS += -f elf64
  else
    ASMFLAGS += -f elf32
  endif
endif
# Flags shared across systems
ifeq ($(ARCH), x86_64)
  ASMFLAGS += -DARCH_X86_64=1
else
  ASMFLAGS += -DARCH_X86_64=0
endif

# Do not use variable length arrays because they don't work in Visual Studio 2013.
WARNINGS = -Wall -Wtype-limits -Wvla

INCLUDEDIRS = -I. -I./strategies -I./extras -I..

STRATEGIES = generic sse2 sse41 avx2 x86_asm altivec

#Detect if cc is gcc for the link time optimization flag
GCCVERSION = $(shell $(CC) --version | grep GCC)
ifneq (, $(GCCVERSION))
  CFLAGS += -flto
  LDFLAGS += -flto -O3
endif
CFLAGS += -std=gnu99 $(INCLUDEDIRS) $(WARNINGS)
LDFLAGS += -lm -pthread
LD = gcc
OBJS = \
  bitstream.o \
  cabac.o \
  checkpoint.o \
  cli.o \
  config.o \
  context.o \
  cu.o \
  encoder.o \
  encoderstate.o \
  rate_control.o \
  filter.o \
  inter.o \
  intra.o \
  kvazaar.o \
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
  strategies/strategies-dct.o \
  strategies/strategies-ipol.o \
  strategies/generic/nal-generic.o \
  strategies/generic/picture-generic.o \
  strategies/sse2/picture-sse2.o \
  strategies/sse41/picture-sse41.o \
  strategies/altivec/picture-altivec.o \
  strategies/avx2/picture-avx2.o \
  strategies/x86_asm/picture-x86-asm.o \
  strategies/generic/dct-generic.o \
  strategies/avx2/dct-avx2.o \
  strategies/generic/ipol-generic.o \
  strategies/avx2/ipol-avx2.o

MAIN_OBJS := encmain.o interface_main.o

TESTDIR := ../tests
TEST_OBJS := \
  $(TESTDIR)/dct_tests.o \
  $(TESTDIR)/intra_sad_tests.o \
  $(TESTDIR)/sad_tests.o \
  $(TESTDIR)/satd_tests.o \
  $(TESTDIR)/speed_tests.o \
  $(TESTDIR)/tests_main.o \
  $(TESTDIR)/test_strategies.o

  
ASMFLAGS += $(INCLUDEDIRS)
ASMOBJS =

AS = yasm

# Compile asm files by default if yasm is present.
ifndef KVZ_DISABLE_ASM
  has_as := $(shell type $(AS) 2>/dev/null)
  ifeq ($(has_as),)
    KVZ_DISABLE_ASM := 1
  endif
endif


PROG  = ./kvazaar
DEBUG = ./kvazaar_debug
TESTS = ./kvazaar_tests

PROGS = $(PROG) $(DEBUG) $(TESTS)

RELEASE_OBJS = $(MAIN_OBJS) $(OBJS) $(ASMOBJS)
DEBUG_OBJS = $(RELEASE_OBJS:.o=_debug.o)
TESTS_OBJS = $(TEST_OBJS) $(OBJS) $(ASMOBJS)


DEPS = $(RELEASE_OBJS:.o=.d) $(DEBUG_OBJS:.o=.d) $(TESTS_OBJS:.o=.d)

GREATEST = ../greatest/greatest.h

all: $(PROG)
.PHONY: all clean build_tests tests debug


# Compile files in strategy directories with appropriate flags.
ifeq ($(ARCH), ppc64)
  strategies/altivec/%.o: EXTRA_FLAGS += -maltivec -fno-lto
else
  strategies/sse2/%.o: EXTRA_FLAGS += -msse2 -fno-lto
  strategies/sse41/%.o: EXTRA_FLAGS += -msse4.1 -fno-lto
  #Needs to be defined on Travis
  ifndef KVZ_DISABLE_AVX2
    strategies/avx2/%.o: EXTRA_FLAGS += -mavx2 -fno-lto
  endif
endif

ifndef KVZ_DISABLE_ASM
  ASMOBJS += strategies/x86_asm/picture-x86-asm-sad.o
  ASMOBJS += strategies/x86_asm/picture-x86-asm-satd.o
  strategies/x86_asm/%.o: EXTRA_FLAGS += -DKVZ_COMPILE_ASM
endif

debug: LDFLAGS := $(filter-out -O3 -O2 -flto, $(LDFLAGS))
debug: CFLAGS := $(filter-out -O3 -O2 -flto, $(CFLAGS))
debug: $(DEBUG)

tests: build_tests
	$(TESTS)

build_tests: CFLAGS := $(filter-out -Werror, $(CFLAGS))
build_tests: init_submodules $(TESTS)

REMOVE_FILES = $(RELEASE_OBJS) $(DEBUG_OBJS) $(TESTS_OBJS) $(DEPS) $(PROGS)
 
$(PROG): $(RELEASE_OBJS)
	$(LD) $^ $(LDFLAGS) -o $@
	
$(DEBUG): $(DEBUG_OBJS)
	$(LD) $^ $(LDFLAGS) -o $@
	
$(TESTS): $(TESTS_OBJS)
	$(LD) $^ $(LDFLAGS) -o $@

%.o %_debug.o: %.asm Makefile
	$(AS) $(ASMFLAGS) -o $@ $<
	
%.o %_debug.o: %.c Makefile
	$(CC) $(CFLAGS) $(EXTRA_FLAGS) -MMD -MP -c $< -o $@

init_submodules:
	git submodule init
	git submodule update

clean:
	rm -f $(REMOVE_FILES)

ifneq ($(MAKECMDGOALS),clean)
-include $(DEPS)
endif
