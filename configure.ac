
AC_INIT([kvazaar], m4_esyscmd([printf $(awk '/#define KVZ_VERSION/ { print $3 }' src/global.h)]))
AC_CONFIG_SRCDIR([src/encmain.c])

# Library version number, modify:
# - When modifying kvazaar.h.
#   - Modify either major or minor.
# - When making a new release.
#   - If major or minor did not change since last release.
#      - Check git history to see if someone forgot to increment major or minor.
#      - Increment release.
#
# major:
#   - Increment when ABI changes, meaning lib users need to be recompiled.
#   - ABI changes when anything existing gets modified, including sizes of structs.
# minor:
#   - Increment when only API changes in a backwards compatible way without breaking ABI.
#   - We count adding parameters to bottom of kvz_config as ABI compatible, because user
#     shouldn't copy that struct or care about it's size.
#   - If not sure, increment major instead.
# release:
#   - Increment when making new releases and major or minor was not changed since last release.
#
# Here is a somewhat sane guide to lib versioning: http://apr.apache.org/versioning.html
ver_major=3
ver_minor=4
ver_release=0

# not used, but it prevents configure from adding a lot of defines to the CFLAGS
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([build-aux])

AM_INIT_AUTOMAKE([-Wall -Werror dist-bzip2 dist-xz foreign subdir-objects])
AM_SILENT_RULES([yes])

AC_PROG_CC
AC_PROG_CC_C99
AM_PROG_AR

# Get fread that can read more than 2GB on 32 bit systems.
AC_SYS_LARGEFILE

LT_INIT([win32-dll])

AX_CHECK_COMPILE_FLAG([-mavx2],   [flag_avx2="true"])
AX_CHECK_COMPILE_FLAG([-msse4.1], [flag_sse4_1="true"])
AX_CHECK_COMPILE_FLAG([-msse2],   [flag_sse2="true"])

AM_CONDITIONAL([HAVE_AVX2], [test x"$flag_avx2" = x"true"])
AM_CONDITIONAL([HAVE_SSE4_1], [test x"$flag_sse4_1" = x"true"])
AM_CONDITIONAL([HAVE_SSE2], [test x"$flag_sse2" = x"true"])

AX_PTHREAD
CFLAGS="-Wall -Wtype-limits -Wvla -I$srcdir/src -ftree-vectorize -fvisibility=hidden $PTHREAD_CFLAGS $CFLAGS"
CPPFLAGS="-DKVZ_DLL_EXPORTS $CPPFLAGS"

AC_SEARCH_LIBS([log], [m c], [], [exit 1])
AC_SEARCH_LIBS([pow], [m c], [], [exit 1])
AC_SEARCH_LIBS([sqrt], [m c], [], [exit 1])
LIBS="$PTHREAD_LIBS $LIBS"


# --disable-werror
AC_ARG_ENABLE([werror], [AS_HELP_STRING([--disable-werror], [don't treat warnings as errors [no]])],
              [], [CFLAGS="-Werror $CFLAGS"]
)


# check for getopt
AC_CHECK_HEADER([getopt.h], [], [CFLAGS="$CFLAGS -I$srcdir/src/extras"])


# host and cpu specific settings
AS_CASE([$host_cpu],
        [i?86], [BITS="32" ASFLAGS="$ASFLAGS -DARCH_X86_64=0" X86="true"],
        [x86_64], [BITS="64" ASFLAGS="$ASFLAGS -DARCH_X86_64=1 -m amd64" X86="true"],
        [powerpc*], [PPC="true"]
)

AS_CASE([$host_os],
        [darwin*], [
         ASFLAGS="$ASFLAGS -f macho$BITS -DPREFIX"
        ],
        [cygwin*|msys*|mingw*], [
         CFLAGS="$CFLAGS -D__USE_MINGW_ANSI_STDIO=1"
         AS_IF(
               [test "x$BITS" = "x32"], [
                ASFLAGS="$ASFLAGS -fwin32 -DPREFIX -DHAVE_ALIGNED_STACK=0"
               ], [
                ASFLAGS="$ASFLAGS -fwin64 -DHAVE_ALIGNED_STACK=1"
               ]
         )
        ],
        [linux*|*kfreebsd*], [
         ASFLAGS="$ASFLAGS -f elf$BITS"
         LDFLAGS="$LDFLAGS -Wl,-z,noexecstack"
         LIBS="$LIBS -lrt"
        ], [
         ASFLAGS="$ASFLAGS -f elf$BITS"
         LDFLAGS="$LDFLAGS -Wl,-z,noexecstack"
        ]
)


# YASM checks
AS_IF([test "x$X86" = "xtrue"], [
       AC_CHECK_TOOL([YASM], [yasm], [no])
])
AS_IF([test "x$YASM" != "xno"], [have_yasm="yes"])

AC_ARG_ENABLE([asm], [AS_HELP_STRING([--disable-asm], [disable assembly [no]])],
              [], [enable_asm="yes"]
)
AS_IF([test "x$enable_asm" != "xno" -a $have_yasm != "yes"],
      [enable_asm="no"]
)


AM_CONDITIONAL([HAVE_X86], [test "x$X86" = "xtrue"])
AM_CONDITIONAL([HAVE_PPC], [test "x$PPC" = "xtrue"])
AM_CONDITIONAL([HAVE_ARM], [test "x$ARM" = "xtrue"])
AM_CONDITIONAL([ENABLE_ASM], [test "x$enable_asm" = "xyes" -a "x$have_yasm" = "xyes" ])

AC_ARG_VAR([ASFLAGS], [ASFLAGS to use for assembler])
AC_SUBST([ASFLAGS])

KVZ_API_VERSION="$ver_major:$ver_minor:$ver_release"
AC_SUBST([KVZ_API_VERSION])


AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/kvazaar.pc
                 tests/Makefile])
AC_OUTPUT

