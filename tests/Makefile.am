
TESTS = $(check_PROGRAMS)

check_PROGRAMS = kvazaar_tests

kvazaar_tests_SOURCES = \
	dct_tests.c \
	intra_sad_tests.c \
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

