# Simple Makefile for Kvazaar HEVC encoder

all: cli lib-shared
.PHONY: all

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

NAME   = libkvazaar
STATIC = $(NAME).a
SO     = $(NAME).so
LIB    = $(SO).$(VER_MAJOR).$(VER_MINOR).$(VER_RELEASE)
DYLIB  = $(NAME).$(VER_MAJOR).dylib
DLL    = kvazaar.dll
IMPLIB = libkvazaar.dll.a
INC    = kvazaar.h \
         kvazaar_version.h

# Compilers and other tools
AS = yasm
LD = $(CC)
INSTALL = install
LN_S = ln -s
RANLIB = ranlib

# Do not use variable length arrays because they don't work in Visual Studio 2013.
WARNINGS = -Werror -Wall -Wtype-limits -Wvla
INCLUDEDIRS = -I. -I./strategies -I./extras -I..

# nasm requires trailing slashes in include directories
ASFLAGS += $(patsubst %,%/,$(INCLUDEDIRS))

# Notes about optimization levels:
# Overall, kvazaar seems to run 15-20 % faster when compiled with gcc.
# gcc-5.1.1: no benefit from -flto, slight harm from -O3
# clang-4.3.0: slight benefit from -O3, more from -flto
CFLAGS += -O2 -g -ftree-vectorize

CFLAGS += -fvisibility=hidden -std=gnu99
CFLAGS += -DKVZ_DLL_EXPORTS
CFLAGS += $(INCLUDEDIRS) $(WARNINGS)
LDFLAGS += -O2 -fvisibility=hidden -lm -pthread

ARCH ?= $(shell uname -m)

ifneq ($(findstring Windows, $(OS)),)
  SYSTEM ?= Windows
else
  SYSTEM ?= $(shell uname -s)
endif

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

INSTALL_TARGETS = install-prog install-pc

# Windows (cygwin/mingw) specific flags
ifeq ($(SYSTEM), Windows)
  ifeq ($(ARCH), x86_64)
    ASFLAGS += -f win64
    ASFLAGS += -DHAVE_ALIGNED_STACK=1
  else
    ASFLAGS += -f win32
    ASFLAGS += -DPREFIX
    ASFLAGS += -DHAVE_ALIGNED_STACK=0
  endif
  CFLAGS += -D__USE_MINGW_ANSI_STDIO=1
  LIBKVAZAAR_SHARED = $(DLL)
  INSTALL_TARGETS += install-dll

# OS X specific flags
else ifeq ($(SYSTEM), Darwin)
  ifeq ($(ARCH), x86_64)
    ASFLAGS += -f macho64
  else
    ASFLAGS += -f macho32
  endif
  ASFLAGS += -DPREFIX
  LIBKVAZAAR_SHARED = $(DYLIB)
  INSTALL_TARGETS += install-dylib

# Default to Linux/elf specific flags
else
  LIBS += -lrt
  LDFLAGS += -Wl,-z,noexecstack
  ifeq ($(ARCH), x86_64)
    ASFLAGS += -f elf64
  else
    ASFLAGS += -f elf32
  endif
  LIBKVAZAAR_SHARED = $(LIB)
  INSTALL_TARGETS += install-lib install-static
endif

# Flags shared across systems
ifeq ($(ARCH), x86_64)
  ASFLAGS += -DARCH_X86_64=1
else
  ASFLAGS += -DARCH_X86_64=0
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
  strategies/altivec/%.o:  EXTRA_FLAGS += -maltivec -fno-lto
  strategies/altivec/%.lo: EXTRA_FLAGS += -maltivec -fno-lto
else
  strategies/sse2/%.o:   EXTRA_FLAGS += -msse2   -fno-lto
  strategies/sse41/%.o:  EXTRA_FLAGS += -msse4.1 -fno-lto
  strategies/sse2/%.lo:  EXTRA_FLAGS += -msse2   -fno-lto
  strategies/sse41/%.lo: EXTRA_FLAGS += -msse4.1 -fno-lto
  # Needs to be defined on Travis
  ifndef KVZ_DISABLE_AVX2
    strategies/avx2/%.o:  EXTRA_FLAGS += -mavx2 -fno-lto
    strategies/avx2/%.lo: EXTRA_FLAGS += -mavx2 -fno-lto
  endif
endif

OBJS = \
  bitstream.o \
  cabac.o \
  checkpoint.o \
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
  search_intra.o \
  search_inter.o \
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

MAIN_OBJS := \
    encmain.o \
    cli.o \
    yuv_io.o

RELEASE_OBJS = $(MAIN_OBJS) $(OBJS)
# Compile separate PIC objects for shared libraries
SHARED_OBJS = $(OBJS:.o=.lo)
DEBUG_OBJS = $(RELEASE_OBJS:.o=_debug.o)
TESTS_OBJS = $(TEST_OBJS) $(OBJS)

DEPS = $(RELEASE_OBJS:.o=.d) $(DEBUG_OBJS:.o=.d) $(TESTS_OBJS:.o=.d)

.PHONY: cli lib-shared lib-static debug tests build_tests
cli: $(PROG)
lib-shared: $(LIBKVAZAAR_SHARED)
lib-static: $(STATIC)

debug: LDFLAGS := $(filter-out -O3 -O2 -flto, $(LDFLAGS))
debug: CFLAGS := $(filter-out -O3 -O2 -flto, $(CFLAGS))
debug: $(DEBUG)

tests: build_tests
	./$(TESTS)

build_tests: CFLAGS := $(filter-out -Werror, $(CFLAGS))
build_tests: init_submodules $(TESTS)

$(LIB):   LDFLAGS += -shared -Wl,-soname,$(SO).$(VER_MAJOR)
$(DLL):   LDFLAGS += -shared -Wl,--out-implib,$(IMPLIB) -o $@
$(DYLIB): LDFLAGS += -dynamiclib \
                     -current_version $(VER_MAJOR).$(VER_MINOR).$(VER_RELEASE) \
                     -compatibility_version $(VER_MAJOR) \
                     -install_name $(LIBDIR)/$@

$(PROG): $(MAIN_OBJS) $(STATIC)
	$(LD) $^ $(LDFLAGS) $(LIBS) -o $@

$(STATIC): $(OBJS)
	$(RM) $@
	$(AR) rc $@ $^
	$(RANLIB) $@

$(LIB) $(DLL) $(DYLIB): $(SHARED_OBJS)
	$(LD) $^ $(LDFLAGS) $(LIBS) -o $@

$(DEBUG): $(DEBUG_OBJS)
	$(LD) $^ $(LDFLAGS) $(LIBS) -o $@

$(TESTS): $(TESTS_OBJS)
	$(LD) $^ $(LDFLAGS) $(LIBS) -o $@


%.o: %.asm Makefile
	$(AS) $(ASFLAGS) -o $@ $<
%_debug.o: %.asm Makefile
	$(AS) $(ASFLAGS) -o $@ $<

# Hopefully the ASM code is PIC. If not, we have to disable for the .so and .dylib.
%.lo: %.asm Makefile
	$(AS) $(ASFLAGS) -o $@ $<
%_debug.lo: %.asm Makefile
	$(AS) $(ASFLAGS) -o $@ $<

# Object files not meant for shared libs need KVZ_STATIC_LIB defined so
# that the names of the API functions don't get mangled on Windows.
%.o: CFLAGS += -DKVZ_STATIC_LIB
%.o: %.c Makefile
	$(CC) $(CFLAGS) $(CPPFLAGS) $(EXTRA_FLAGS) -MMD -MP -c $< -o $@
%_debug.o: %.c Makefile
	$(CC) $(CFLAGS) $(CPPFLAGS) $(EXTRA_FLAGS) -MMD -MP -c $< -o $@

ifneq ($(SYSTEM), Windows)
  # Only use -fpic when not on Windows
  %.lo: CFLAGS += -fpic
endif
%.lo: %.c Makefile
	$(CC) $(CFLAGS) $(CPPFLAGS) $(EXTRA_FLAGS) -MMD -MP -c $< -o $@
%_debug.lo: %.c Makefile
	$(CC) $(CFLAGS) $(CPPFLAGS) $(EXTRA_FLAGS) -MMD -MP -c $< -o $@


.PHONY: kvazaar.pc init_submodules install install-pc install-prog install-lib
.PHONY: install-dylib install-dll clean
kvazaar.pc: kvazaar.pc.in Makefile
	sed -e "s;@prefix@;$(PREFIX);" -e "s;@libdir@;$(LIBDIR);" \
            -e "s;@VERSION@;$(VER_MAJOR).$(VER_MINOR).$(VER_RELEASE);" \
            -e "s;@LIBS@;$(LIBS);" kvazaar.pc.in > $@

init_submodules:
	git submodule init
	git submodule update

install: $(INSTALL_TARGETS)

install-pc: kvazaar.pc
	$(INSTALL) -d $(DESTDIR)$(LIBDIR)/pkgconfig
	$(INSTALL) -m644 $^ $(DESTDIR)$(LIBDIR)/pkgconfig

install-prog: $(PROG)
	$(INSTALL) -d $(DESTDIR)$(BINDIR)
	$(INSTALL) -m755 $(PROG) $(DESTDIR)$(BINDIR)

install-static: $(STATIC)
	$(INSTALL) -d $(DESTDIR)$(INCDIR)
	$(INSTALL) -d $(DESTDIR)$(LIBDIR)
	$(INSTALL) -m644 $(INC) $(DESTDIR)$(INCDIR)
	$(INSTALL) -m644 $(STATIC) $(DESTDIR)$(LIBDIR)

install-lib: $(LIB)
	$(INSTALL) -d $(DESTDIR)$(INCDIR)
	$(INSTALL) -d $(DESTDIR)$(LIBDIR)
	$(INSTALL) -m644 $(INC) $(DESTDIR)$(INCDIR)
	$(INSTALL) -m644 $(LIB) $(DESTDIR)$(LIBDIR)
	$(LN_S) -f -T $(LIB) $(DESTDIR)$(LIBDIR)/$(SO).$(VER_MAJOR)
	$(LN_S) -f -T $(LIB) $(DESTDIR)$(LIBDIR)/$(SO)

install-dylib: $(DYLIB)
	$(INSTALL) -d $(DESTDIR)$(INCDIR)
	$(INSTALL) -d $(DESTDIR)$(LIBDIR)
	$(INSTALL) -m644 $(INC) $(DESTDIR)$(INCDIR)
	$(INSTALL) -m644 $(DYLIB) $(DESTDIR)$(LIBDIR)
	$(LN_S) -f $(DYLIB) $(DESTDIR)$(LIBDIR)/$(NAME).dylib

install-dll: $(DLL)
	$(INSTALL) -d $(DESTDIR)$(DLLDIR)
	$(INSTALL) -d $(DESTDIR)$(INCDIR)
	$(INSTALL) -d $(DESTDIR)$(LIBDIR)
	$(INSTALL) -m644 $(DLL) $(DESTDIR)$(DLLDIR)
	$(INSTALL) -m644 $(INC) $(DESTDIR)$(INCDIR)
	$(INSTALL) -m644 $(IMPLIB) $(DESTDIR)$(LIBDIR)

clean:
	$(RM) $(RELEASE_OBJS) $(SHARED_OBJS) $(DEBUG_OBJS) $(TESTS_OBJS) $(DEPS)
	$(RM) $(PROG) $(DEBUG) $(TESTS) $(STATIC) $(LIB) $(DLL) $(IMPLIB) $(DYLIB)

ifneq ($(MAKECMDGOALS),clean)
-include $(DEPS)
endif

# Disable default rules
.SUFFIXES:

.PHONY: all clean build_tests tests debug init_submodules
.PHONY: install install-prog install-lib install-dll
