ACLOCAL_AMFLAGS = -I m4


bin_PROGRAMS = kvazaar
lib_LTLIBRARIES = libkvazaar.la


dist_doc_DATA = COPYING CREDITS README.md

EXTRA_DIST = \
	build \
	doc \
	docs.doxy \
	greatest \
	SConstruct \
	tests \
	tools \
	src/Makefile \
	src/SConscript \
	src/extras/getopt.c \
	src/extras/getopt.h


pkgconfigdir = $(libdir)/pkgconfig
nodist_pkgconfig_DATA = src/kvazaar.pc


include_HEADERS = \
	src/kvazaar.h \
	src/kvazaar_version.h

noinst_HEADERS = \
	src/extras/x86inc.asm


kvazaar_SOURCES = \
	src/encmain.c \
	src/cli.h \
	src/cli.c \
	src/yuv_io.c \
	src/yuv_io.h
kvazaar_LDADD = libkvazaar.la $(LIBS)


libkvazaar_la_SOURCES = \
	src/bitstream.c \
	src/bitstream.h \
	src/cabac.c \
	src/cabac.h \
	src/checkpoint.c \
	src/checkpoint.h \
	src/config.c \
	src/config.h \
	src/context.c \
	src/context.h \
	src/cu.c \
	src/cu.h \
	src/encoder.c \
	src/encoder.h \
	src/encoderstate.c \
	src/encoderstate.h \
	src/encoder_state-bitstream.c \
	src/encoder_state-bitstream.h \
	src/encoder_state-ctors_dtors.c \
	src/encoder_state-ctors_dtors.h \
	src/encoder_state-geometry.c \
	src/encoder_state-geometry.h \
	src/filter.c \
	src/filter.h \
	src/global.h \
	src/image.c \
	src/image.h \
	src/imagelist.c \
	src/imagelist.h \
	src/input_frame_buffer.c \
	src/input_frame_buffer.h \
	src/inter.c \
	src/interface_main.c \
	src/inter.h \
	src/intra.c \
	src/intra.h \
	src/kvazaar.c \
	src/kvazaar_internal.h \
	src/nal.c \
	src/nal.h \
	src/rate_control.c \
	src/rate_control.h \
	src/rdo.c \
	src/rdo.h \
	src/sao.c \
	src/sao.h \
	src/scalinglist.c \
	src/scalinglist.h \
	src/search.c \
	src/search.h \
	src/search_inter.c \
	src/search_inter.h \
	src/search_intra.c \
	src/search_intra.h \
	src/tables.c \
	src/tables.h \
	src/threadqueue.c \
	src/threadqueue.h \
	src/threads.h \
	src/transform.c \
	src/transform.h \
	src/videoframe.c \
	src/videoframe.h \
	src/strategies/altivec/picture-altivec.c \
	src/strategies/altivec/picture-altivec.h \
	src/strategies/generic/dct-generic.c \
	src/strategies/generic/dct-generic.h \
	src/strategies/generic/intra-generic.c \
	src/strategies/generic/intra-generic.h \
	src/strategies/generic/ipol-generic.c \
	src/strategies/generic/ipol-generic.h \
	src/strategies/generic/nal-generic.c \
	src/strategies/generic/nal-generic.h \
	src/strategies/generic/picture-generic.c \
	src/strategies/generic/picture-generic.h \
	src/strategies/generic/quant-generic.c \
	src/strategies/generic/quant-generic.h \
	src/strategies/sse2/picture-sse2.c \
	src/strategies/sse2/picture-sse2.h \
	src/strategies/sse41/picture-sse41.c \
	src/strategies/sse41/picture-sse41.h \
	src/strategies/strategies-common.h \
	src/strategies/strategies-dct.c \
	src/strategies/strategies-dct.h \
	src/strategies/strategies-intra.c \
	src/strategies/strategies-intra.h \
	src/strategies/strategies-ipol.c \
	src/strategies/strategies-ipol.h \
	src/strategies/strategies-nal.c \
	src/strategies/strategies-nal.h \
	src/strategies/strategies-picture.c \
	src/strategies/strategies-picture.h \
	src/strategies/strategies-quant.c \
	src/strategies/strategies-quant.h \
	src/strategies/x86_asm/picture-x86-asm.c \
	src/strategies/x86_asm/picture-x86-asm.h \
	src/strategyselector.c \
	src/strategyselector.h

avx2_sources = \
	src/strategies/avx2/dct-avx2.c \
	src/strategies/avx2/dct-avx2.h \
	src/strategies/avx2/intra-avx2.c \
	src/strategies/avx2/intra-avx2.h \
	src/strategies/avx2/ipol-avx2.c \
	src/strategies/avx2/ipol-avx2.h \
	src/strategies/avx2/picture-avx2.c \
	src/strategies/avx2/picture-avx2.h \
	src/strategies/avx2/quant-avx2.c \
	src/strategies/avx2/quant-avx2.h

libkvazaar_la_SOURCES += $(avx2_sources)
libkvazaar_la_LDFLAGS = $(AM_LDFLAGS) -no-undefined -version-number $(KVZ_API_VERSION)


if HAVE_PPC

src/strategies/altivec/picture-altivec.lo: CFLAGS += -maltivec

else #HAVE_PPC

if HAVE_X86
src/strategies/sse2/picture-sse2.lo: CFLAGS += -msse2
src/strategies/sse41/picture-sse41.lo: CFLAGS += -msse4.1
$(avx2_sources): CFLAGS += -mavx2

if ENABLE_ASM
libkvazaar_la_SOURCES += \
	src/strategies/x86_asm/picture-x86-asm-sad.asm \
	src/strategies/x86_asm/picture-x86-asm-sad.h \
	src/strategies/x86_asm/picture-x86-asm-satd.asm \
	src/strategies/x86_asm/picture-x86-asm-satd.h

src/strategies/x86_asm/picture-x86-asm.lo: CFLAGS += -DKVZ_COMPILE_ASM
src/strategies/x86_asm/picture-x86-asm-sad.lo: src/strategies/x86_asm/picture-x86-asm-sad.asm
src/strategies/x86_asm/picture-x86-asm-satd.lo: src/strategies/x86_asm/picture-x86-asm-satd.asm
endif #ENABLE_ASM
endif #HAVE_X86

endif #HAVE_PPC


yasm_verbose = $(yasm_verbose_@AM_V@)
yasm_verbose_ = $(yasm_verbose_@AM_DEFAULT_V@)
yasm_verbose_0 = @echo "  YASM    " $@;

.asm.lo:
	$(yasm_verbose)$(LIBTOOL) --mode=compile $(YASM) $(ASFLAGS) $< -o $@ -prefer-non-pic 1>/dev/null

