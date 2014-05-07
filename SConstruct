""" SConstruct file fore building Kvazaar with Scons.

This file defines two targets, x86 and x64, and builds construction
environments for both.

TODO:
- add debug builds
- whole program optimization for gcc

"""

import os
import platform

Help("""
Type: 'scons x86' to build a 32-bit release version,
     'scons x64' to build a 64-bit release version.
     'scons x86 x64 -c' to clear all build directories.
""")


# Create construction environments for 32 and 64 bit builds.
# Visual studio needs the architecture to be set in this stage. It can not be
# modified later.
env_x86 = Environment(
        #tools=['mingw'],
        ASCOM='yasm $ASFLAGS -o $TARGET $SOURCES',
        ENV={'PATH': os.environ['PATH']},  # to find yasm on Windows
        TARGET_ARCH='x86',  # for Visual Studio
        )
env_x64 = Environment(
        #tools=['mingw'],
        ASCOM='yasm $ASFLAGS -o $TARGET $SOURCES',
        ENV={'PATH': os.environ['PATH']}, # to find yasm on Windows
        TARGET_ARCH='amd64',  # for Visual Studio
        )

# YASM flags for different architectures. The object file format and name
# mangling must be the same as used by the C compiler for that OS.
# Indexed by platform.system().
yasm_flags = {
    'Windows': {
        'x86': '-f win32 -DPREFIX',
        'x64': '-f win64'},
    'Darwin': {
        'x86': '-f macho32 -DPREFIX',
        'x64': '-f macho64 -DPREFIX'},
    'Linux': {  # Flags for Unix-like.
        'x86': '-f elf32',
        'x64': '-f elf64'},
    'all': {  # Flags for all systems.
        'x86': ' -DARCH_X86_64=0 -m x86',
        'x64': ' -DARCH_X86_64=1 -m amd64'},
}
# Set yasm flags for OS and architecture.
target_yasm_flags = yasm_flags.get(platform.system(), yasm_flags['Linux'])
env_x86.Replace(ASFLAGS=target_yasm_flags['x86'])
env_x64.Replace(ASFLAGS=target_yasm_flags['x64'])
env_x86.Append(ASFLAGS=yasm_flags['all']['x86'])
env_x64.Append(ASFLAGS=yasm_flags['all']['x64'])


# Try and deal with all the remaining compiler specific differences that I
# really with scons would handle.
if 'MSVS' in env_x86:
    # /MD = multithreaded DLL runtime
    # /Ox = full optimization, sets /Ob2, /Og, /Oi, /Ot, /Oy
    # /GL = enable whole program optimization
    # /LTCG = link time code generation
    # /arch:SSE2 = use SSE2 (x86 only)
    env_x86.Append(
            CCFLAGS='/MD /Ox /GL /arch:SSE2',
            LINKFLAGS='/LTCG')
    env_x64.Append(
            CCFLAGS='/MD /Ox /GL /openmp',
            LINKFLAGS='/LTCG')
else:
    # GCC flags
    # -m for arch, -O2 for optimization, -lm for math lib
    env_x86.MergeFlags('-m32 -O2 -lm -march=native -fopenmp')
    env_x64.MergeFlags('-m64 -O2 -lm -march=native -fopenmp')

# VS2010 linker and mingw64 need TMP.
if 'TMP' in os.environ:
    env_x86['ENV']['TMP'] = os.environ['TMP']
    env_x64['ENV']['TMP'] = os.environ['TMP']

env_x86.MergeFlags('-Iextras')
env_x64.MergeFlags('-Iextras')
    
preprocessor_defines = ARGUMENTS.get('D', '')
if preprocessor_defines:
    for define in preprocessor_defines.split():
        env_x86.MergeFlags('-D' + define)
        env_x64.MergeFlags('-D' + define)

# Declare build targets.
x86 = SConscript('src/SConscript',
        exports={'env': env_x86},
        variant_dir='scons_build_x86',
        duplicate=False)
Alias('x86', x86)

x64 = SConscript('src/SConscript',
        exports={'env': env_x64},
        variant_dir='scons_build_x64',
        duplicate=False)
Alias('x64', x64)


if platform.machine().endswith('64'):
    Default(x64)
else:
    Default(x86)

