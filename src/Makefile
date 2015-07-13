# Simple Makefile for Kvazaar HEVC encoder

# Installation locations
DESTDIR =
PREFIX  = /usr/local
BINDIR  = $(PREFIX)/bin
INCDIR  = $(PREFIX)/include
LIBDIR  = $(PREFIX)/lib
DLLDIR  = $(BINDIR)

# Library version number
VER_MAJOR = 0
VER_MINOR = 0
VER_RELEASE = 0

PROG   = kvazaar
DEBUG  = kvazaar_debug
TESTS  = kvazaar_tests
SO     = libkvazaar.so
LIB    = $(SO).$(VER_MAJOR).$(VER_MINOR).$(VER_RELEASE)
DLL    = kvazaar.dll
IMPLIB = libkvazaar.dll.a
INC    = kvazaar.h \
         kvazaar_version.h

# Compilers and other tools
AS = yasm
LD = gcc
INSTALL = install
LN_S = ln -s

# Do not use variable length arrays because they don't work in Visual Studio 2013.
WARNINGS = -Wall -Wtype-limits -Wvla
INCLUDEDIRS = -I. -I./strategies -I./extras -I..

ASFLAGS += $(INCLUDEDIRS)
CFLAGS += -O2 -g -Werror -ftree-vectorize -fpic -std=gnu99
CFLAGS += $(INCLUDEDIRS) $(WARNINGS)
LDFLAGS += -lm -pthread

ifeq (, $(ARCH))
  ARCH = $(shell uname -m)
endif
SYSTEM = $(shell uname -s)

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
    ASFLAGS += -f win64
    ASFLAGS += -DHAVE_ALIGNED_STACK=1
  else
    ASFLAGS += -f win32
    ASFLAGS += -DPREFIX
    ASFLAGS += -DHAVE_ALIGNED_STACK=0
  endif
  CFLAGS += -D__USE_MINGW_ANSI_STDIO=1
# OS X specific flags
else ifeq ($(SYSTEM),Darwin)
  ifeq ($(ARCH), x86_64)
    ASFLAGS += -f macho64
  else
    ASFLAGS += -f macho32
  endif
  ASFLAGS += -DPREFIX
# Default to Linux/elf specific flags
else
  LDFLAGS += -lrt
  ifeq ($(ARCH), x86_64)
    ASFLAGS += -f elf64
  else
    ASFLAGS += -f elf32
  endif
endif
# Flags shared across systems
ifeq ($(ARCH), x86_64)
  ASFLAGS += -DARCH_X86_64=1
else
  ASFLAGS += -DARCH_X86_64=0
endif

#Detect if cc is gcc for the link time optimization flag
GCCVERSION = $(shell $(CC) --version | grep GCC)
ifneq (, $(GCCVERSION))
  CFLAGS += -flto
  LDFLAGS += -flto -O3
endif

# Compile asm files by default if yasm is present.
ifndef KVZ_DISABLE_ASM
  has_as := $(shell type $(AS) 2>/dev/null)
  ifeq ($(has_as),)
    KVZ_DISABLE_ASM := 1
  endif
endif

ifdef CHECKPOINTS
  CFLAGS += -DCHECKPOINTS
endif

# Compile files in strategy directories with appropriate flags.
ifeq ($(ARCH), ppc64)
  strategies/altivec/%.o: EXTRA_FLAGS += -maltivec -fno-lto
else
  strategies/sse2/%.o: EXTRA_FLAGS += -msse2 -fno-lto
  strategies/sse41/%.o: EXTRA_FLAGS += -msse4.1 -fno-lto
  # Needs to be defined on Travis
  ifndef KVZ_DISABLE_AVX2
    strategies/avx2/%.o: EXTRA_FLAGS += -mavx2 -fno-lto
  endif
endif

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
  yuv_io.o \
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

ifndef KVZ_DISABLE_ASM
  strategies/x86_asm/%.o: EXTRA_FLAGS += -DKVZ_COMPILE_ASM

  OBJS += \
    strategies/x86_asm/picture-x86-asm-sad.o \
    strategies/x86_asm/picture-x86-asm-satd.o
endif

TESTDIR := ../tests
TEST_OBJS := \
  $(TESTDIR)/dct_tests.o \
  $(TESTDIR)/intra_sad_tests.o \
  $(TESTDIR)/sad_tests.o \
  $(TESTDIR)/satd_tests.o \
  $(TESTDIR)/speed_tests.o \
  $(TESTDIR)/tests_main.o \
  $(TESTDIR)/test_strategies.o

MAIN_OBJS := encmain.o interface_main.o

RELEASE_OBJS = $(MAIN_OBJS) $(OBJS)
DEBUG_OBJS = $(RELEASE_OBJS:.o=_debug.o)
TESTS_OBJS = $(TEST_OBJS) $(OBJS)

DEPS = $(RELEASE_OBJS:.o=.d) $(DEBUG_OBJS:.o=.d) $(TESTS_OBJS:.o=.d)

ifneq ( ,$(findstring Windows, $(OS)))
  all: $(PROG) $(DLL)
else
  all: $(PROG) $(LIB)
endif

debug: LDFLAGS := $(filter-out -O3 -O2 -flto, $(LDFLAGS))
debug: CFLAGS := $(filter-out -O3 -O2 -flto, $(CFLAGS))
debug: $(DEBUG)

tests: build_tests
	$(TESTS)

build_tests: CFLAGS := $(filter-out -Werror, $(CFLAGS))
build_tests: init_submodules $(TESTS)

$(PROG): $(RELEASE_OBJS)
	$(LD) $^ $(LDFLAGS) -o $@

$(LIB): $(RELEASE_OBJS)
	$(LD) $^ $(LDFLAGS) -shared -Wl,-soname,$(SO).$(VER_MAJOR) -o $@

$(DLL): $(RELEASE_OBJS)
	$(LD) $^ $(LDFLAGS) -shared -Wl,--out-implib,$(IMPLIB) -o $@

$(DEBUG): $(DEBUG_OBJS)
	$(LD) $^ $(LDFLAGS) -o $@

$(TESTS): $(TESTS_OBJS)
	$(LD) $^ $(LDFLAGS) -o $@

%.o %_debug.o: %.asm Makefile
	$(AS) $(ASFLAGS) -o $@ $<

%.o %_debug.o: %.c Makefile
	$(CC) $(CFLAGS) $(EXTRA_FLAGS) -MMD -MP -c $< -o $@

init_submodules:
	git submodule init
	git submodule update

ifneq ( ,$(findstring Windows, $(OS)))
  install: install-prog install-dll
else
  install: install-prog install-lib
endif

install-prog: $(PROG)
	$(INSTALL) -d $(DESTDIR)$(BINDIR)
	$(INSTALL) -m755 $(PROG) -t $(DESTDIR)$(BINDIR)

install-lib: $(LIB)
	$(INSTALL) -d $(DESTDIR)$(INCDIR)
	$(INSTALL) -d $(DESTDIR)$(LIBDIR)
	$(INSTALL) -m644 $(INC) -t $(DESTDIR)$(INCDIR)
	$(INSTALL) -m644 $(LIB) -t $(DESTDIR)$(LIBDIR)
	$(LN_S) -f -T $(LIB) $(DESTDIR)$(LIBDIR)/$(SO).$(VER_MAJOR)
	$(LN_S) -f -T $(LIB) $(DESTDIR)$(LIBDIR)/$(SO)

install-dll: $(DLL)
	$(INSTALL) -d $(DESTDIR)$(DLLDIR)
	$(INSTALL) -d $(DESTDIR)$(INCDIR)
	$(INSTALL) -d $(DESTDIR)$(LIBDIR)
	$(INSTALL) -m644 $(DLL) -t $(DESTDIR)$(DLLDIR)
	$(INSTALL) -m644 $(INC) -t $(DESTDIR)$(INCDIR)
	$(INSTALL) -m644 $(IMPLIB) -t $(DESTDIR)$(LIBDIR)

clean:
	$(RM) $(RELEASE_OBJS) $(DEBUG_OBJS) $(TESTS_OBJS) $(DEPS)
	$(RM) $(PROG) $(DEBUG) $(TESTS) $(LIB) $(DLL) $(IMPLIB)

ifneq ($(MAKECMDGOALS),clean)
-include $(DEPS)
endif

# Disable default rules
.SUFFIXES:

.PHONY: all clean build_tests tests debug init_submodules
.PHONY: install install-prog install-lib install-dll
