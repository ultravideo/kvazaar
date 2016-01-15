
bin_PROGRAMS = kvazaar
lib_LTLIBRARIES = libkvazaar.la


EXTRA_DIST = \
	extras/getopt.c \
	extras/getopt.h


pkgconfigdir = $(libdir)/pkgconfig
nodist_pkgconfig_DATA = kvazaar.pc


include_HEADERS = \
	kvazaar.h \
	kvazaar_version.h

noinst_HEADERS = \
	extras/x86inc.asm


kvazaar_SOURCES = \
	encmain.c \
	cli.h \
	cli.c \
	yuv_io.c \
	yuv_io.h
kvazaar_LDADD = libkvazaar.la $(LIBS)


libkvazaar_la_SOURCES = \
	bitstream.c \
	bitstream.h \
	cabac.c \
	cabac.h \
	checkpoint.c \
	checkpoint.h \
	config.c \
	config.h \
	context.c \
	context.h \
	cu.c \
	cu.h \
	encoder.c \
	encoder.h \
	encoderstate.c \
	encoderstate.h \
	encoder_state-bitstream.c \
	encoder_state-bitstream.h \
	encoder_state-ctors_dtors.c \
	encoder_state-ctors_dtors.h \
	encoder_state-geometry.c \
	encoder_state-geometry.h \
	filter.c \
	filter.h \
	global.h \
	image.c \
	image.h \
	imagelist.c \
	imagelist.h \
	input_frame_buffer.c \
	input_frame_buffer.h \
	inter.c \
	inter.h \
	intra.c \
	intra.h \
	kvazaar.c \
	kvazaar_internal.h \
	nal.c \
	nal.h \
	rate_control.c \
	rate_control.h \
	rdo.c \
	rdo.h \
	sao.c \
	sao.h \
	scalinglist.c \
	scalinglist.h \
	search.c \
	search.h \
	search_inter.c \
	search_inter.h \
	search_intra.c \
	search_intra.h \
	tables.c \
	tables.h \
	threadqueue.c \
	threadqueue.h \
	threads.h \
	transform.c \
	transform.h \
	videoframe.c \
	videoframe.h \
	strategies/altivec/picture-altivec.c \
	strategies/altivec/picture-altivec.h \
	strategies/generic/dct-generic.c \
	strategies/generic/dct-generic.h \
	strategies/generic/intra-generic.c \
	strategies/generic/intra-generic.h \
	strategies/generic/ipol-generic.c \
	strategies/generic/ipol-generic.h \
	strategies/generic/nal-generic.c \
	strategies/generic/nal-generic.h \
	strategies/generic/picture-generic.c \
	strategies/generic/picture-generic.h \
	strategies/generic/quant-generic.c \
	strategies/generic/quant-generic.h \
	strategies/sse2/picture-sse2.c \
	strategies/sse2/picture-sse2.h \
	strategies/sse41/picture-sse41.c \
	strategies/sse41/picture-sse41.h \
	strategies/strategies-common.h \
	strategies/strategies-dct.c \
	strategies/strategies-dct.h \
	strategies/strategies-intra.c \
	strategies/strategies-intra.h \
	strategies/strategies-ipol.c \
	strategies/strategies-ipol.h \
	strategies/strategies-nal.c \
	strategies/strategies-nal.h \
	strategies/strategies-picture.c \
	strategies/strategies-picture.h \
	strategies/strategies-quant.c \
	strategies/strategies-quant.h \
	strategies/x86_asm/picture-x86-asm.c \
	strategies/x86_asm/picture-x86-asm.h \
	strategyselector.c \
	strategyselector.h

avx2_sources = \
	strategies/avx2/dct-avx2.c \
	strategies/avx2/dct-avx2.h \
	strategies/avx2/intra-avx2.c \
	strategies/avx2/intra-avx2.h \
	strategies/avx2/ipol-avx2.c \
	strategies/avx2/ipol-avx2.h \
	strategies/avx2/picture-avx2.c \
	strategies/avx2/picture-avx2.h \
	strategies/avx2/quant-avx2.c \
	strategies/avx2/quant-avx2.h

libkvazaar_la_SOURCES += $(avx2_sources)
libkvazaar_la_LDFLAGS = $(AM_LDFLAGS) -no-undefined -version-number $(KVZ_API_VERSION)


if HAVE_PPC

strategies/altivec/picture-altivec.lo: CFLAGS += -maltivec

else #HAVE_PPC

if HAVE_X86
strategies/sse2/picture-sse2.lo: CFLAGS += -msse2
strategies/sse41/picture-sse41.lo: CFLAGS += -msse4.1
$(avx2_sources): CFLAGS += -mavx2

if ENABLE_ASM
libkvazaar_la_SOURCES += \
	strategies/x86_asm/picture-x86-asm-sad.asm \
	strategies/x86_asm/picture-x86-asm-sad.h \
	strategies/x86_asm/picture-x86-asm-satd.asm \
	strategies/x86_asm/picture-x86-asm-satd.h

strategies/x86_asm/picture-x86-asm.lo: CFLAGS += -DKVZ_COMPILE_ASM
strategies/x86_asm/picture-x86-asm-sad.lo: strategies/x86_asm/picture-x86-asm-sad.asm
strategies/x86_asm/picture-x86-asm-satd.lo: strategies/x86_asm/picture-x86-asm-satd.asm
endif #ENABLE_ASM
endif #HAVE_X86

endif #HAVE_PPC


yasm_verbose = $(yasm_verbose_@AM_V@)
yasm_verbose_ = $(yasm_verbose_@AM_DEFAULT_V@)
yasm_verbose_0 = @echo "  YASM    " $@;

.asm.lo:
	$(yasm_verbose)$(LIBTOOL) --mode=compile $(YASM) -I$(srcdir)/extras $(ASFLAGS) $< -o $@ -prefer-non-pic 1>/dev/null

