""" SConstruct file fore building Kvazaar with Scons.

This is an additional cross platform way to build the program.
The main way is to use the Makefile and that is always kept up to date.

This mostly exists so it can be used with a custom SConscript file that
builds any old version of Kvazaar, ensuring that the compilation settings are
the same and any compilation from a specific version of Kvazaar will not
affect the comparison. The SConscript included in git version only compiles
the current version.

"""

import os
import platform


vars = Variables()

default_arch = 'amd64' if platform.machine().endswith('64') else 'x86'
vars.Add(EnumVariable('arch', 'Set target arch.', default_arch,
        allowed_values=('x86', 'x64', 'amd64', 'ppc'),
        map={'x64': 'amd64'}))

vars.Add(PathVariable('win32pthreads',
        'Path to win32-pthreads dir.',
        r'./../pthreads.2'))

vars.Add(PathVariable('copyto',
        'Copy exe and required DLLs to this dir.',
        '',
        PathVariable.PathAccept))

vars.Add(BoolVariable('dump_env',
        'Dump of construction environment to stdout.',
        False))

vars.Add(BoolVariable('use_yasm',
        'Use yasm.',
        True))


# Visual studio needs the architecture to be set in this stage. It can not be
# modified later. Other tools ignore TARGET_ARCH.
# The variable substitution is done in the environment construction so $arch
# will not work here. Get the value manually.
arch = ARGUMENTS.get('arch', default_arch)
if arch == 'x64':
  arch = 'amd64'
env = Environment(
        variables=vars,
        tools=['msvc', 'mslink', 'nasm'],
        ENV={'PATH': os.environ['PATH']},  # to find yasm
        TARGET_ARCH=arch,  # for Visual Studio
        AS='yasm',
        )


Help("""
Example: 'scons arch=x64' to compile for amd64.
         'scons --jobs=8' to compile in parallel.
""" + vars.GenerateHelpText(env))

if 'MSVS' in env:
    compiler = 'msvs'
elif 'MSYSTEM' in os.environ:
	compiler = 'mingw'
else:
	compiler = 'gcc'

env['use_yasm'] = env['use_yasm'] and env['arch'] not in ('ppc',)

# pthreads_dll location for copyto
pthreads_dll = None

# Try and deal with all the remaining compiler specific differences that I
# really wish scons would handle.
if compiler in ('msvs',):
    arch_dir = {'x86': r'x86', 'amd64': r'x64'}[env['arch']]
    pthreads_dll = os.path.join(env['win32pthreads'], 'dll', arch_dir, 'pthreadVC2.dll')
    
    env.Append(
            CCFLAGS=r'/MD /Ox /GL /wd"4028"',
            LINKFLAGS=r'/LTCG',
            LIBPATH=r'#$win32pthreads\lib\\' + arch_dir,
            LIBS=r'pthreadVC2',
            CPPPATH=r'#$win32pthreads\include')
else:
    env.MergeFlags('-O2 -pthread -lm -lrt -march=native')
    if env['arch'] == 'x86':
        env.MergeFlags('-m32')
    elif env['arch'] == 'amd64':
        env.MergeFlags('-m64')
    
    # platform.system() lies on msys2, so just try to detect msys/mingw.
    if 'MSYSTEM' in os.environ:
        # __USE_MINGW_ANSI_STDIO required on mingw for printf to function.
        env.MergeFlags('-D' + '__USE_MINGW_ANSI_STDIO=1')


# VS2010 linker and mingw64 need TMP.
if compiler in ('msvs', 'mingw'):
	if 'TMP' in os.environ:
	    env['ENV']['TMP'] = os.environ['TMP']

env.MergeFlags('-I. -Iextras -Istrategies')

# Take comma separated list from 'D=' and pass them on to
# preprocessor as defines.
preprocessor_defines = ARGUMENTS.get('D', '')
if preprocessor_defines:
    for define in preprocessor_defines.split(','):
        env.MergeFlags('-D' + define)

# Declare build targets.
variant_dir = 'scons_build'


class EnvinronmentContainer(object):
    """Class for initializing and holding optimization environments.
    
    As far as I know, making a separate environment is the only way to compile
    different objects with different flags in scons. So this constructs all
    the required environments.
    
    Yasm is also its own environment although it's not strictly necessary.
    It does however allow KVZ_COMPILE_ASM to only those files that require
    it, avoiding a complete recompile.
    
    """
    
    def __init__(self, env, compiler, arch):
        # If optimization is not supported for arch, use default env.
        self.env = env
        self.sse2 = env
        self.sse41 = env
        self.avx = env
        self.avx2 = env
        self.altivec = env
        self.asm = env
        
        if compiler == 'msvs':
            self.avx = env.Clone()
            self.avx.Append(CCFLAGS='/arch:AVX')
            self.avx2 = env.Clone()
            self.avx2.Append(CCFLAGS='/arch:AVX2')
        elif compiler in ('gcc',) and arch in ('x86', 'x64'):
            self.sse2 = env.Clone()
            self.sse2.Append(CCFLAGS='-msse2')
            self.sse41 = env.Clone().Append(CCFLAGS='-msse4.1')
            self.avx = env.Clone().Append(CCFLAGS='-mavx')
            self.avx2 = env.Clone().Append(CCFLAGS='-mavx2')
        elif compiler in ('gcc',) and arch in ('ppc'):
            self.altivec = env.Clone().Append(CCFLAGS='-maltivec')

        if env['use_yasm']:
            self._init_yasm()
    
    
    def _init_yasm(self):
        self.asm = env.Clone()
        self.asm.MergeFlags('-D' + 'KVZ_COMPILE_ASM')
        
        # YASM flags for different architectures. The object file format and name
        # mangling must be the same as used by the C compiler for that OS.
        # Indexed by platform.system().
        yasm_flags = {
            'Windows': {
                'x86':   '-f win32 -DPREFIX -DHAVE_ALIGNED_STACK=0 ',
                'amd64': '-f win64 -DHAVE_ALIGNED_STACK=1',
            },
            'Darwin': {
                'x86':   '-f macho32 -DPREFIX',
                'amd64': '-f macho64 -DPREFIX',
            },
            'Linux': {  # Flags for Unix-like.
                'x86':   '-f elf32',
                'amd64': '-f elf64',
            },
            'all': {  # Flags for all systems.
                'x86':   ' -I./src/extras -DARCH_X86_64=0 -m x86',
                'amd64': ' -I./src/extras -DARCH_X86_64=1 -m amd64',
            },
        }

        # Set yasm flags for OS and architecture.
        target_yasm_flags = yasm_flags.get(platform.system(), yasm_flags['Linux'])
        self.asm.Replace(ASFLAGS=target_yasm_flags[env['arch']])
        self.asm.Append(ASFLAGS=yasm_flags['all'][env['arch']])
        self.asm.Replace(AS='yasm')


envs = EnvinronmentContainer(env, compiler, env['arch'])


program = SConscript('src/SConscript',
        exports={'envs': envs},
        variant_dir=variant_dir,
        duplicate=False)


# Copy exe and dll's to the path from 'copyto=' argument.
copy_dir = env['copyto']
if copy_dir:
    env.Install(copy_dir, program)
    if pthreads_dll:
        env.Install(copy_dir, pthreads_dll)
    env.Alias('copyto', copy_dir)
    Default([program, 'copyto'])
else:
    Default(program)


# Dumpo environment to stdout.
if env['dump_env']:
    import difflib
    
    def get_diff(orig_env, new_env):
        return "\n".join(difflib.unified_diff(orig_env, new_env.Dump().splitlines(), n=0))
    
    env_dump = env.Dump().splitlines()
    print "== env"
    print "\n".join(env_dump)
    print "== diff env envs.sse2"
    print get_diff(env_dump, envs.sse2)
    print "== diff env envs.sse41"
    print get_diff(env_dump, envs.sse41)
    print "== diff env envs.avx"
    print get_diff(env_dump, envs.avx)
    print "== diff env envs.avx2"
    print get_diff(env_dump, envs.avx2)
    print "== diff env envs.altivec"
    print get_diff(env_dump, envs.altivec)
    print "== diff env envs.asm"
    print get_diff(env_dump, envs.asm)
    
