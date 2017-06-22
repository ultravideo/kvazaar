
TESTS = $(check_PROGRAMS) \
    test_external_symbols.bash \
    test_gop.bash \
    test_interlace.bash \
    test_intra.bash \
    test_invalid_input.bash \
    test_mv_constraint.bash \
    test_owf_wpp_tiles.bash \
    test_rate_control.bash \
    test_slices.bash \
    test_smp.bash \
    test_tools.bash \
    test_weird_shapes.bash

check_PROGRAMS = kvazaar_tests

kvazaar_tests_SOURCES = \
	dct_tests.c \
	intra_sad_tests.c \
	mv_cand_tests.c \
	sad_tests.c \
	sad_tests.h \
	satd_tests.c \
	satd_tests.h \
	speed_tests.c \
	tests_main.c \
	test_strategies.c \
	test_strategies.h
kvazaar_tests_CFLAGS = -I$(srcdir) -I$(top_srcdir) -I$(top_srcdir)/src
kvazaar_tests_LDFLAGS = -static $(top_builddir)/src/libkvazaar.la $(LIBS)
