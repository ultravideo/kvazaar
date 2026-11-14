bin_PROGRAMS = kvazaar
lib_LTLIBRARIES = libkvazaar.la


EXTRA_DIST = \
	extras/getopt.c \
	extras/getopt.h \
	version.h.in \
	kvazaarCMake.pc.in \
	threadwrapper

pkgconfigdir = $(libdir)/pkgconfig
nodist_pkgconfig_DATA = kvazaar.pc


include_HEADERS = \
	kvazaar.h

noinst_LTLIBRARIES = \
	libaltivec.la \
	libavx2.la \
	libneon.la \
	libsse2.la \
	libsse41.la

kvazaar_SOURCES = \
	encmain.c \
	cli.h \
	cli.c \
	yuv_io.c \
	yuv_io.h

kvazaar_LDADD = libkvazaar.la $(LIBS)

kvazaar_CPPFLAGS = -DKVZ_VERSION="`$(srcdir)/../tools/version.sh`"

if USE_CRYPTOPP
kvazaar_LINK = $(LIBTOOL) $(AM_V_lt) --tag=CXX $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=link $(CXXLD) $(AM_CXXFLAGS) $(CXXFLAGS) \
	$(AM_LDFLAGS) $(LDFLAGS) -o $@
else
kvazaar_LINK = $(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=link $(CCLD) $(AM_CFLAGS) $(CFLAGS) \
	$(AM_LDFLAGS) $(LDFLAGS) -o $@
endif

# Sources of the library that are compiled for every bit depth. The SIMD
# sources are in separate libraries because they need their own flags.
KVZ_LIB_SOURCES = \
	multilib_rename.h \
	bitstream.c \
	bitstream.h \
	cabac.c \
	cabac.h \
	checkpoint.c \
	checkpoint.h \
	cfg.c \
	cfg.h \
	constraint.c \
	constraint.h \
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
	encode_coding_tree.c \
	encode_coding_tree.h \
	fast_coeff_cost.c \
	fast_coeff_cost.h \
	filter.c \
	filter.h \
	global.h \
	gop.h \
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
	kvz_math.h \
	ml_intra_cu_depth_pred.c \
	ml_intra_cu_depth_pred.h \
	nal.c \
	nal.h \
	rate_control.c \
	rate_control.h \
	reconstruct.c \
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
	sei.h \
	tables.c \
	tables.h \
	threadqueue.c \
	threadqueue.h \
	threads.h \
	transform.c \
	transform.h \
	videoframe.c \
	videoframe.h \
	work_tree_pool.c \
	work_tree_pool.h \
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
	strategies/generic/sao-generic.c \
	strategies/generic/sao-generic.h \
	strategies/generic/encode_coding_tree-generic.c \
	strategies/generic/encode_coding_tree-generic.h \
	strategies/missing-intel-intrinsics.h \
	strategies/optimized_sad_func_ptr_t.h \
	strategies/generic/sao_shared_generics.h \
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
	strategies/strategies-sao.c \
	strategies/strategies-sao.h \
	strategies/strategies-encode.c \
	strategies/strategies-encode.h \
	strategyselector.c \
	strategyselector.h \
	extras/libmd5.c \
	extras/libmd5.h \
	extras/crypto.h

KVZ_ALTIVEC_SOURCES = \
	strategies/altivec/picture-altivec.c \
	strategies/altivec/picture-altivec.h

KVZ_AVX2_SOURCES = \
	strategies/strategies-cpp.h \
	strategies/avx2/avx2_common_functions.h \
	strategies/avx2/dct-avx2.cpp \
	strategies/avx2/dct-avx2.h \
	strategies/avx2/intra-avx2.cpp \
	strategies/avx2/intra-avx2.h \
	strategies/avx2/ipol-avx2.c \
	strategies/avx2/ipol-avx2.h \
	strategies/avx2/picture-avx2.cpp \
	strategies/avx2/picture-avx2.h \
	strategies/avx2/quant-avx2.c \
	strategies/avx2/quant-avx2.h \
	strategies/avx2/reg_sad_pow2_widths-avx2.h \
	strategies/avx2/sao-avx2.c \
	strategies/avx2/sao-avx2.h \
	strategies/avx2/encode_coding_tree-avx2.c \
	strategies/avx2/encode_coding_tree-avx2.h

KVZ_NEON_SOURCES = \
	strategies/strategies-cpp.h \
	strategies/neon/neon-common.h \
	strategies/neon/dct-neon.cpp \
	strategies/neon/dct-neon.h \
	strategies/neon/intra-neon.cpp \
	strategies/neon/intra-neon.h \
	strategies/neon/ipol-neon.cpp \
	strategies/neon/ipol-neon.h \
	strategies/neon/picture-neon.cpp \
	strategies/neon/picture-neon.h \
	strategies/neon/quant-neon.c \
	strategies/neon/quant-neon.h \
	strategies/neon/sao-neon.c \
	strategies/neon/sao-neon.h

KVZ_SSE2_SOURCES = \
	strategies/sse2/picture-sse2.c \
	strategies/sse2/picture-sse2.h

KVZ_SSE41_SOURCES = \
	strategies/sse41/picture-sse41.c \
	strategies/sse41/picture-sse41.h \
	strategies/sse41/reg_sad_pow2_widths-sse41.h

# Flags for SIMD strategies written in C++ (see strategies/strategies-cpp.h).
# They use templates only: no exceptions, RTTI or C++ runtime, so the library
# is still linked as C. Symbols are hidden like in the C code.
KVZ_STRATEGY_CXXFLAGS = $(AM_CXXFLAGS) -std=c++17 -fno-exceptions -fno-rtti -fvisibility=hidden \
	-Wall -Wextra -Wvla -Wno-sign-compare -Wno-unused-parameter -ftree-vectorize \
	-I$(srcdir) -I$(srcdir)/extras

if HAVE_PPC
if HAVE_ALTIVEC
KVZ_ALTIVEC_CFLAGS = -maltivec
endif
endif #HAVE_PPC

if HAVE_X86
if HAVE_AVX2_GCC
KVZ_AVX2_CFLAGS = -mavx2 -mbmi -mabm -mpopcnt -mlzcnt -mbmi2
endif
if HAVE_AVX2_CLANG
KVZ_AVX2_CFLAGS = -mavx2 -mbmi -mpopcnt -mlzcnt -mbmi2
endif
if HAVE_SSE4_1
KVZ_SSE41_CFLAGS = -msse4.1
endif
if HAVE_SSE2
KVZ_SSE2_CFLAGS = -msse2
endif
endif #HAVE_X86

# The default bit depth (the lowest one given to configure --with-bit-depths)
# is compiled as usual and its kvz_api_get returns the APIs of the others
# (KVZ_MULTILIB_<n>BIT).
KVZ_DEFAULT_CPPFLAGS = -DKVZ_BIT_DEPTH=$(KVZ_DEFAULT_BIT_DEPTH)

libkvazaar_la_SOURCES = $(KVZ_LIB_SOURCES)
libkvazaar_la_CPPFLAGS = $(KVZ_DEFAULT_CPPFLAGS) $(KVZ_MULTILIB_CPPFLAGS)

libkvazaar_la_LIBADD = \
	libaltivec.la \
	libavx2.la \
	libneon.la \
	libsse2.la \
	libsse41.la

libkvazaar_la_LDFLAGS = $(AM_LDFLAGS) -no-undefined -version-number $(KVZ_API_VERSION)

if USE_CRYPTOPP
libkvazaar_la_SOURCES += \
	extras/crypto.h \
	extras/crypto.cpp
libkvazaar_la_LINK = $(LIBTOOL) $(AM_V_lt) --tag=CXX $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=link $(CXXLD) $(AM_CXXFLAGS) $(CXXFLAGS) \
	$(libkvazaar_la_LDFLAGS) $(LDFLAGS) -o $@
else
libkvazaar_la_LINK = $(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=link $(CCLD) $(AM_CFLAGS) $(CFLAGS) \
	$(libkvazaar_la_LDFLAGS) $(LDFLAGS) -o $@
endif

libaltivec_la_SOURCES = $(KVZ_ALTIVEC_SOURCES)
libaltivec_la_CPPFLAGS = $(KVZ_DEFAULT_CPPFLAGS)
libaltivec_la_CFLAGS = $(KVZ_ALTIVEC_CFLAGS)

libavx2_la_SOURCES = $(KVZ_AVX2_SOURCES)
libavx2_la_CPPFLAGS = $(KVZ_DEFAULT_CPPFLAGS)
libavx2_la_CFLAGS = $(KVZ_AVX2_CFLAGS)
libavx2_la_CXXFLAGS = $(KVZ_STRATEGY_CXXFLAGS) $(KVZ_AVX2_CFLAGS)

# The AArch64 NEON strategies need no extra compiler flags: NEON is part of
# the AArch64 baseline. On other architectures they compile to nothing.
libneon_la_SOURCES = $(KVZ_NEON_SOURCES)
libneon_la_CPPFLAGS = $(KVZ_DEFAULT_CPPFLAGS)
libneon_la_CXXFLAGS = $(KVZ_STRATEGY_CXXFLAGS)

libsse2_la_SOURCES = $(KVZ_SSE2_SOURCES)
libsse2_la_CPPFLAGS = $(KVZ_DEFAULT_CPPFLAGS)
libsse2_la_CFLAGS = $(KVZ_SSE2_CFLAGS)

libsse41_la_SOURCES = $(KVZ_SSE41_SOURCES)
libsse41_la_CPPFLAGS = $(KVZ_DEFAULT_CPPFLAGS)
libsse41_la_CFLAGS = $(KVZ_SSE41_CFLAGS)

# Each additional bit depth is a copy of the library compiled with
# KVZ_BIT_DEPTH set and multilib_rename.h included in every file, which
# renames the external symbols (e.g. kvz_api_get -> kvz_api_get_10bit).

KVZ_10BIT_CPPFLAGS = -DKVZ_BIT_DEPTH=10 -include $(srcdir)/multilib_rename.h

if KVZ_MULTILIB_10BIT
noinst_LTLIBRARIES += \
	libkvazaar_10bit.la \
	libaltivec_10bit.la \
	libavx2_10bit.la \
	libneon_10bit.la \
	libsse2_10bit.la \
	libsse41_10bit.la
libkvazaar_la_LIBADD += libkvazaar_10bit.la
endif

libkvazaar_10bit_la_SOURCES = $(KVZ_LIB_SOURCES)
libkvazaar_10bit_la_CPPFLAGS = $(KVZ_10BIT_CPPFLAGS)
libkvazaar_10bit_la_LIBADD = \
	libaltivec_10bit.la \
	libavx2_10bit.la \
	libneon_10bit.la \
	libsse2_10bit.la \
	libsse41_10bit.la

libaltivec_10bit_la_SOURCES = $(KVZ_ALTIVEC_SOURCES)
libaltivec_10bit_la_CPPFLAGS = $(KVZ_10BIT_CPPFLAGS)
libaltivec_10bit_la_CFLAGS = $(KVZ_ALTIVEC_CFLAGS)

libavx2_10bit_la_SOURCES = $(KVZ_AVX2_SOURCES)
libavx2_10bit_la_CPPFLAGS = $(KVZ_10BIT_CPPFLAGS)
libavx2_10bit_la_CFLAGS = $(KVZ_AVX2_CFLAGS)
libavx2_10bit_la_CXXFLAGS = $(KVZ_STRATEGY_CXXFLAGS) $(KVZ_AVX2_CFLAGS)

libneon_10bit_la_SOURCES = $(KVZ_NEON_SOURCES)
libneon_10bit_la_CPPFLAGS = $(KVZ_10BIT_CPPFLAGS)
libneon_10bit_la_CXXFLAGS = $(KVZ_STRATEGY_CXXFLAGS)

libsse2_10bit_la_SOURCES = $(KVZ_SSE2_SOURCES)
libsse2_10bit_la_CPPFLAGS = $(KVZ_10BIT_CPPFLAGS)
libsse2_10bit_la_CFLAGS = $(KVZ_SSE2_CFLAGS)

libsse41_10bit_la_SOURCES = $(KVZ_SSE41_SOURCES)
libsse41_10bit_la_CPPFLAGS = $(KVZ_10BIT_CPPFLAGS)
libsse41_10bit_la_CFLAGS = $(KVZ_SSE41_CFLAGS)

KVZ_12BIT_CPPFLAGS = -DKVZ_BIT_DEPTH=12 -include $(srcdir)/multilib_rename.h

if KVZ_MULTILIB_12BIT
noinst_LTLIBRARIES += \
	libkvazaar_12bit.la \
	libaltivec_12bit.la \
	libavx2_12bit.la \
	libneon_12bit.la \
	libsse2_12bit.la \
	libsse41_12bit.la
libkvazaar_la_LIBADD += libkvazaar_12bit.la
endif

libkvazaar_12bit_la_SOURCES = $(KVZ_LIB_SOURCES)
libkvazaar_12bit_la_CPPFLAGS = $(KVZ_12BIT_CPPFLAGS)
libkvazaar_12bit_la_LIBADD = \
	libaltivec_12bit.la \
	libavx2_12bit.la \
	libneon_12bit.la \
	libsse2_12bit.la \
	libsse41_12bit.la

libaltivec_12bit_la_SOURCES = $(KVZ_ALTIVEC_SOURCES)
libaltivec_12bit_la_CPPFLAGS = $(KVZ_12BIT_CPPFLAGS)
libaltivec_12bit_la_CFLAGS = $(KVZ_ALTIVEC_CFLAGS)

libavx2_12bit_la_SOURCES = $(KVZ_AVX2_SOURCES)
libavx2_12bit_la_CPPFLAGS = $(KVZ_12BIT_CPPFLAGS)
libavx2_12bit_la_CFLAGS = $(KVZ_AVX2_CFLAGS)
libavx2_12bit_la_CXXFLAGS = $(KVZ_STRATEGY_CXXFLAGS) $(KVZ_AVX2_CFLAGS)

libneon_12bit_la_SOURCES = $(KVZ_NEON_SOURCES)
libneon_12bit_la_CPPFLAGS = $(KVZ_12BIT_CPPFLAGS)
libneon_12bit_la_CXXFLAGS = $(KVZ_STRATEGY_CXXFLAGS)

libsse2_12bit_la_SOURCES = $(KVZ_SSE2_SOURCES)
libsse2_12bit_la_CPPFLAGS = $(KVZ_12BIT_CPPFLAGS)
libsse2_12bit_la_CFLAGS = $(KVZ_SSE2_CFLAGS)

libsse41_12bit_la_SOURCES = $(KVZ_SSE41_SOURCES)
libsse41_12bit_la_CPPFLAGS = $(KVZ_12BIT_CPPFLAGS)
libsse41_12bit_la_CFLAGS = $(KVZ_SSE41_CFLAGS)
