Kvazaar
=======
An open-source HEVC encoder licensed under LGPLv2.1

Join channel #kvazaar_hevc in Freenode IRC network to contact us.

Kvazaar is not yet finished and does not implement all the features of
HEVC. Compression performance will increase as we add more coding tools.

http://ultravideo.cs.tut.fi/#encoder for more information.

[![Build Status](https://travis-ci.org/ultravideo/kvazaar.svg?branch=master)](https://travis-ci.org/ultravideo/kvazaar)

## Using Kvazaar

[comment]: # (BEGIN KVAZAAR HELP MESSAGE)
```
Usage:
kvazaar -i <input> --input-res <width>x<height> -o <output>

Required:
  -i, --input                : Input file
      --input-res <res>      : Input resolution [auto]
                               auto: detect from file name
                               <int>x<int>: width times height
  -o, --output               : Output file

Presets:
      --preset=<preset>      : Set options to a preset [medium]
                                   - ultrafast, superfast, veryfast, faster,
                                     fast, medium, slow, slower, veryslow
                                     placebo

Input:
  -n, --frames <integer>     : Number of frames to code [all]
      --seek <integer>       : First frame to code [0]
      --input-fps <num>/<denom> : Framerate of the input video [25.0]
      --source-scan-type <string> : Set source scan type [progressive].
                                   - progressive: progressive scan
                                   - tff: top field first
                                   - bff: bottom field first
      --input-format         : P420 or P400
      --input-bitdepth       : 8-16
      --loop-input           : Re-read input file forever

Options:
      --help                 : Print this help message and exit
      --version              : Print version information and exit
      --aud                  : Use access unit delimiters
      --debug <string>       : Output encoders reconstruction.
      --cpuid <integer>      : Disable runtime cpu optimizations with value 0.
      --hash                 : Decoded picture hash [checksum]
                                   - none: 0 bytes
                                   - checksum: 18 bytes
                                   - md5: 56 bytes
      --no-psnr              : Don't calculate PSNR for frames
      --no-info              : Don't add encoder info SEI.

Video structure:
  -q, --qp <integer>         : Quantization Parameter [32]
  -p, --period <integer>     : Period of intra pictures [0]
                               - 0: only first picture is intra
                               - 1: all pictures are intra
                               - 2-N: every Nth picture is intra
      --vps-period <integer> : Specify how often the video parameter set is
                               re-sent. [0]
                                   - 0: only send VPS with the first frame
                                   - N: send VPS with every Nth intra frame
  -r, --ref <integer>        : Reference frames, range 1..15 [3]
      --gop <string>         : Definition of GOP structure [0]
                                   - 0: disabled
                                   - 8: B-frame pyramid of length 8
                                   - lp-<string>: lp-gop definition
                                         (e.g. lp-g8d4t2, see README)
      --cqmfile <string>     : Custom Quantization Matrices from a file
      --bitrate <integer>    : Target bitrate. [0]
                                   - 0: disable rate-control
                                   - N: target N bits per second
      --lossless             : Use lossless coding
      --mv-constraint        : Constrain movement vectors
                                   - none: no constraint
                                   - frametile: constrain within the tile
                                   - frametilemargin: constrain even more

Compression tools:
      --deblock [<beta:tc>]  : Deblocking
                                     - beta: between -6 and 6
                                     - tc: between -6 and 6
      --(no-)sao             : Sample Adaptive Offset
      --(no-)rdoq            : Rate-Distortion Optimized Quantization
      --(no-)signhide        : Sign Hiding
      --(no-)smp             : Symmetric Motion Partition
      --(no-)amp             : Asymmetric Motion Partition
      --rd <integer>         : Intra mode search complexity
                                   - 0: skip intra if inter is good enough
                                   - 1: rough intra mode search with SATD
                                   - 2: refine intra mode search with SSE
      --(no-)mv-rdo          : Rate-Distortion Optimized motion vector costs
      --(no-)full-intra-search
                             : Try all intra modes during rough search.
      --(no-)transform-skip  : Transform skip
      --me <string>          : Integer motion estimation
                                   - hexbs: Hexagon Based Search
                                   - tz:    Test Zone Search
                                   - full:  Full Search
                                   - full8, full16, full32, full64
      --subme <integer>      : Set fractional pixel motion estimation level
                                   - 0: only integer motion estimation
                                   - 1: + 1/2-pixel horizontal and vertical
                                   - 2: + 1/2-pixel diagonal
                                   - 3: + 1/4-pixel horizontal and vertical
                                   - 4: + 1/4-pixel diagonal
      --pu-depth-inter <int>-<int>
                             : Range for sizes for inter predictions
                                   - 0, 1, 2, 3: from 64x64 to 8x8
      --pu-depth-intra <int>-<int> : Range for sizes for intra predictions
                                   - 0, 1, 2, 3, 4: from 64x64 to 4x4
      --(no-)bipred          : Bi-prediction
      --(no-)cu-split-termination
                             : CU split search termination condition
                                   - off: Never terminate cu-split search
                                   - zero: Terminate with zero residual
      --(no-)me-early-termination : ME early termination condition
                                   - off: Don't terminate early
                                   - on: Terminate early
                                   - sensitive: Terminate even earlier
      --(no-)implicit-rdpcm  : Implicit residual DPCM
                               Currently only supported with lossless coding.
      --(no-)tmvp            : Temporal Motion Vector Prediction
      --(no-)rdoq-skip       : Skips RDOQ for 4x4 blocks

Parallel processing:
      --threads <integer>    : Number of threads to use [auto]
                                   - 0: process everything with main thread
                                   - N: use N threads for encoding
                                   - auto: select based on number of cores
      --owf <integer>        : Frame parallelism [auto]
                                   - N: Process N-1 frames at a time
                                   - auto: Select automatically
      --(no-)wpp             : Wavefront parallel processing [enabled]
                               Enabling tiles automatically disables WPP.
                               To enable WPP with tiles, re-enable it after
                               enabling tiles.
      --tiles <int>x<int>    : Split picture into width x height uniform tiles.
      --tiles-width-split <string>|u<int> :
                               Specifies a comma separated list of pixel
                               positions of tiles columns separation coordinates.
                               Can also be u followed by and a single int n,
                               in which case it produces columns of uniform width.
      --tiles-height-split <string>|u<int> :
                               Specifies a comma separated list of pixel
                               positions of tiles rows separation coordinates.
                               Can also be u followed by and a single int n,
                               in which case it produces rows of uniform height.

Video Usability Information:
      --sar <width:height>   : Specify Sample Aspect Ratio
      --overscan <string>    : Specify crop overscan setting [undef]
                                   - undef, show, crop
      --videoformat <string> : Specify video format [undef]
                                   - component, pal, ntsc, secam, mac, undef
      --range <string>       : Specify color range [tv]
                                   - tv, pc
      --colorprim <string>   : Specify color primaries [undef]
                                   - undef, bt709, bt470m, bt470bg,
                                     smpte170m, smpte240m, film, bt2020
      --transfer <string>    : Specify transfer characteristics [undef]
                                   - undef, bt709, bt470m, bt470bg,
                                     smpte170m, smpte240m, linear, log100,
                                     log316, iec61966-2-4, bt1361e,
                                     iec61966-2-1, bt2020-10, bt2020-12
      --colormatrix <string> : Specify color matrix setting [undef]
                                   - undef, bt709, fcc, bt470bg, smpte170m,
                                     smpte240m, GBR, YCgCo, bt2020nc, bt2020c
      --chromaloc <integer>  : Specify chroma sample location (0 to 5) [0]

Deprecated parameters: (might be removed at some point)
  -w, --width                 : Use --input-res
  -h, --height                : Use --input-res
```
[comment]: # (END KVAZAAR HELP MESSAGE)

### For example:

    kvazaar -i BQMall_832x480_60.yuv --input-res 832x480 -o out.hevc -n 600 -q 32

The only accepted input format so far is 8-bit YUV 4:2:0.


### LP-GOP syntax
The LP-GOP syntax is "lp-g(num)d(num)r(num)t(num)", where
- g = GOP length.
- d = Number of GOP layers.
- r = Number of references, where one reference is always the previous picture,
  unless temporal scaling is used. The others are key-frames.
- t = How many references to skip for temporal scaling, where 4 means only
  every fourth picture needs to be decoded.



## Presets
The names of the presets are the same as with x264: ultrafast,
superfast, veryfast, faster, fast, medium, slow, slower, veryslow and
placebo. The effects of the presets are listed in the following table,
where the names have been abbreviated to fit the layout in GitHub.

                     | 0-uf  | 1-sf  | 2-vf  | 3-fr  | 4-f   | 5-m   | 6-s   | 7-sr  | 8-vs  | 9-p
-------------------- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | -----
rd                   | 0     | 0     | 0     | 1     | 1     | 1     | 1     | 1     | 1     | 1
pu-depth-intra       | 2-3   | 2-3   | 2-3   | 2-3   | 2-3   | 1-3   | 1-3   | 1-3   | 1-4   | 1-4
pu-depth-inter       | 2-3   | 2-3   | 2-3   | 1-3   | 1-3   | 1-3   | 1-3   | 0-3   | 0-3   | 0-3
me                   | hexbs | hexbs | hexbs | hexbs | hexbs | hexbs | hexbs | hexbs | hexbs | tz
ref                  | 1     | 1     | 1     | 1     | 1     | 1     | 2     | 2     | 3     | 4
deblock              | 1     | 1     | 1     | 1     | 1     | 1     | 1     | 1     | 1     | 1
signhide             | 0     | 0     | 0     | 0     | 0     | 0     | 1     | 1     | 1     | 1
subme                | 0     | 0     | 2     | 2     | 4     | 4     | 4     | 4     | 4     | 4
sao                  | 0     | 1     | 1     | 1     | 1     | 1     | 1     | 1     | 1     | 1
rdoq                 | 0     | 0     | 0     | 0     | 0     | 1     | 1     | 1     | 1     | 1
rdoq-skip            | 1     | 1     | 1     | 1     | 1     | 1     | 1     | 1     | 1     | 0
transform-skip       | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 1
mv-rdo               | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 1
full-intra-search    | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0
smp                  | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 1
amp                  | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 0     | 1
cu-split-termination | zero  | zero  | zero  | zero  | zero  | zero  | zero  | zero  | zero  | off
me-early-termination | sens. | sens. | sens. | sens. | on    | on    | on    | on    | on    | off


## Kvazaar library

See [kvazaar.h](src/kvazaar.h) for the library API and its
documentation.

When using the static Kvazaar library on Windows, macro `KVZ_STATIC_LIB`
must be defined. On other platforms it's not strictly required.

The needed linker and compiler flags can be obtained with pkg-config.

## Compiling Kvazaar

If you have trouble regarding compiling the source code, please make an
[issue](https://github.com/ultravideo/kvazaar/issues) about in Github.
Others might encounter the same problem and there is probably much to
improve in the build process. We want to make this as simple as
possible.

### Required libraries
- For Visual Studio, the pthreads-w32 library is required. Platforms
  with native POSIX thread support don't need anything.
  - The project file expects the library to be in ../pthreads.2/
    relative to Kvazaar. You can just extract the pre-built library
    there.
  - The executable needs pthreadVC2.dll to be present. Either install it
    somewhere or ship it with the executable.

### Autotools

Depending on the platform, some additional tools are required for compiling Kvazaar with autotools.
For Ubuntu, the required packages are `automake autoconf libtool m4 build-essential yasm`.

Run the following commands to compile and install Kvazaar.

    ./autogen.sh
    ./configure
    make
    sudo make install

See `./configure --help` for more options.

### OS X
- The program should compile and work on OS X but you might need a newer
  version of GCC than what comes with the platform.

### Visual Studio
- VS2010 and older do not have support for some of the C99 features that
  we use. Please use VS2013 or newer or GCC (MinGW) to compile on
  Windows.
- Project files can be found under build/.
- Requires external [vsyasm.exe](http://yasm.tortall.net/Download.html)
  in %PATH%
  - Run `rundll32 sysdm.cpl,EditEnvironmentVariables` and add PATH to
    user variables
- Building the Kvazaar library is not yet supported.

### Docker

This project includes a [Dockerfile](./Dockerfile), which enables building for Docker. Kvazaar is also available in the Docker Hub [`ultravideo/kvazaar`](https://hub.docker.com/r/ultravideo/kvazaar/)
Build using Docker: `docker build -t kvazaar .`
Example usage: `docker run -i -a STDIN -a STDOUT kvazaar -i - --input-res=320x240 -o - < testfile_320x240.yuv > out.265`
For other examples, see [Dockerfile](./Dockerfile)

### Visualization (Windows only)

Branch `visualizer` has a visual studio project, which can be compiled to enable visualization feature in Kvazaar.

Additional Requirements: [`SDL2`](https://www.libsdl.org/download-2.0.php), [`SDL2-ttf`](https://www.libsdl.org/projects/SDL_ttf/).

Directory `visualizer_extras` is expected to be found from the same directory level as the kvazaar project directory. Inside should be directories `include` and `lib` found from the development library zip packages.

`SDL2.dll`, `SDL2_ttf.dll`, `libfreetype-6.dll`, `zlib1.dll`, and `pthreadVC2.dll` should be placed in the working directory (i.e. the folder the `kvazaar.exe` is in after compiling the `kvazaar_cli` project/solution) when running the visualizer. The required `.dll` can be found in the aforementioned `lib`-folder (`lib\x64`) and the dll folder inside the pthreads folder (see `Required libraries`).

Note: The solution should be compiled on the x64 platform in visual studio.

Optional font file `arial.ttf` is to be placed in the working directory, if block info tool is used.

## Contributing to Kvazaar

See http://github.com/ultravideo/kvazaar/wiki/List-of-suggested-topics
for a list of topics you might want to examine if you would like to do
something bigger than a bug fix but don't know what yet.


### Code documentation

You can generate Doxygen documentation pages by running the command
"doxygen docs.doxy". Here is a rough sketch of the module structure:
![Kvazaar module hierarchy](https://github.com/ultravideo/kvazaar/blob/master/doc/kvazaar_module_hierarchy.png)

### For version control we try to follow these conventions:

- Master branch always produces a working bitstream (can be decoded with
  HM).
- Commits for new features and major changes/fixes put to a sensibly
  named feature branch first and later merged to the master branch.
- Always merge the feature branch to the master branch, not the other
  way around, with fast-forwarding disabled if necessary. We have found
  that this differentiates between working and unfinished versions
  nicely.
- Every commit should at least compile. Producing a working bitstream is
  nice as well, but not always possible. Features may be temporarily
  disabled to produce a working bitstream, but remember to re-enable
  them before merging to master.


### Testing

- We do not have a proper testing framework yet. We test mainly by
  decoding the bitstream with HM and checking that the result matches
  the encoders own reconstruction.
- You should at least test that HM decodes a bitstream file made with
  your changes without throwing checksum errors. If your changes
  shouldn't alter the bitstream, you should check that they don't.
- We would like to have a suite of automatic tests that also check for
  BD-rate increase and speed decrease in addition to checking that the
  bitstream is valid. As of yet there is no such suite.


### Unit tests
- There are some unit tests located in the tests directory. We would
  like to have more.
- The Visual Studio project links the unit tests against the actual .lib
  file used by the encoder. There is no Makefile as of yet.
- The unit tests use "greatest" unit testing framework. It is included
  as a submodule, but getting it requires the following commands to be
  run in the root directory of kvazaar:

        git submodule init
        git submodule update


### Code style

We try to follow the following conventions:
- C99 without features not supported by Visual Studio 2013 (VLAs).
 - // comments allowed and encouraged.
- Follow overall conventions already established in the code.
- Indent by 2 spaces. (no tabs)
- { on the same line for control logic and on the next line for functions
- Reference and deference next to the variable name.
- Variable names in lowered characters with words divided by underscore.
- Maximum line length 79 characters when possible.
- Functions only used inside the module shouldn't be defined in the
  module header. They can be defined in the beginning of the .c file if
  necessary.


### Resources for HEVC bitstream features

- A good first resource for HEVC bitstream is JCTVC-N1002 High
  Efficiency Video Coding (HEVC) Test Model 12 (HM12) Encoder
  Description
- Many good articles regarding specific parts of HEVC can be found on
  IEEE Transactions on Circuits and Systems for Video Technology,
  Combined issue on High Efficiency Video Coding (HEVC) Standards and
  Research
- The specification tends to follow the reference implementation, not
  the other way around, so check HM if the specification is unclear.
