
AC_INIT([kvazaar], m4_esyscmd([printf $(awk '/#define KVZ_VERSION/ { print $3 }' src/global.h)]))
AC_CONFIG_SRCDIR([src/encmain.c])

# Library version number
ver_major=2
ver_minor=4
ver_release=0

# not used, but it prevents configure from adding a lot of defines to the CFLAGS
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([build-aux])

AM_INIT_AUTOMAKE([-Wall -Werror dist-bzip2 dist-xz foreign subdir-objects])
AM_SILENT_RULES([yes])

AC_PROG_CC
AM_PROG_AR

LT_INIT([win32-dll])


CFLAGS="-Wall -Wtype-limits -Wvla -I$srcdir/src -ftree-vectorize -fvisibility=hidden -std=gnu99 -DKVZ_DLL_EXPORTS $CFLAGS"
LIBS="-lm $LIBS"


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
        [powerpc*], [PPC="true"],
        [arm*], [ARM="true"]
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
        ], [
         ASFLAGS="$ASFLAGS -f elf$BITS"
         CFLAGS="$CFLAGS -pthread"
         LDFLAGS="$LDFLAGS -pthread -Wl,-z,noexecstack"
         LIBS="$LIBS -lrt -lpthread"
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

# use KVZ_API_VERSION for sonames
#KVZ_API_VERSION="$($AWK '/#define KVZ_API_VERSION/ { print $3 }' $srcdir/src/kvazaar_version.h)"

# use library version for sonames (specified in src/Makefile)
KVZ_API_VERSION="$ver_major:$ver_minor:$ver_release"
AC_SUBST([KVZ_API_VERSION])


AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/kvazaar.pc
                 tests/Makefile])
AC_OUTPUT
