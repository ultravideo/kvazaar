#Kvazaar
An open-source HEVC encoder licensed under GPLv2

Join channel #kvazaar_hevc in Freenode IRC network to contact us.

Kvazaar is not yet finished and does not implement all the features of HEVC. Compression performance will increase as we add more coding tools.

http://ultravideo.cs.tut.fi/#encoder for more information.

##Using Kvazaar
Currently most of the features are turned on/off from the code on compile time, but they are
meant to be user configurable later.

    Usage:
    kvazaar -i <input> -w <width> -h <height> -o <output>

    Optional parameters:
          -n, --frames <integer>     : number of frames to code [all]
          -q, --qp <integer>         : Quantization Parameter [32]
          -p, --period <integer>     : Period of intra pictures [0]
                                         0: only first picture is intra
                                         1: all pictures are intra
                                         2-N: every Nth picture is intra
              --no-deblock           : Disable deblocking filter
              --deblock <beta:tc>    : Deblocking filter parameters
                                       beta and tc range is -6..6 [0:0]
              --no-sao               : Disable sample adaptive offset
              --aud                  : Use access unit delimiters
              --cqmfile <string>     : Custom Quantization Matrices from a file

      Video Usability Information:
              --sar <width:height>   : Specify Sample Aspect Ratio
              --overscan <string>    : Specify crop overscan setting ["undef"]
                                         - undef, show, crop
              --videoformat <string> : Specify video format ["undef"]
                                         - component, pal, ntsc, secam, mac, undef
              --range <string>       : Specify color range ["tv"]
                                         - tv, pc
              --colorprim <string>   : Specify color primaries ["undef"]
                                         - undef, bt709, bt470m, bt470bg,
                                           smpte170m, smpte240m, film, bt2020
              --transfer <string>    : Specify transfer characteristics ["undef"]
                                         - undef, bt709, bt470m, bt470bg,
                                           smpte170m, smpte240m, linear, log100,
                                           log316, iec61966-2-4, bt1361e,
                                           iec61966-2-1, bt2020-10, bt2020-12
              --colormatrix <string> : Specify color matrix setting ["undef"]
                                         - undef, bt709, fcc, bt470bg, smpte170m,
                                           smpte240m, GBR, YCgCo, bt2020nc, bt2020c
              --chromaloc <integer>  : Specify chroma sample location (0 to 5) [0]

Example:

    kvazaar -i <INPUT_YUV> -w <WIDTH> -h <HEIGHT> -o <OUTPUT.BIN> -n <NUMBER_OF_FRAMES> -q <QP>

eg. `kvazaar -i BQMall_832x480_60.yuv -w 832 -h 480 -o out.bin -n 600 -q 32`

The only accepted input format so far is 8-bit YUV 4:2:0.


##Compiling Kvazaar

If you have trouble regarding compiling the source code, please make an [issue](https://github.com/ultravideo/kvazaar/issues) about in Github. Others might encounter the same problem and there is probably much to improve in the build process. We want to make this as simple as possible.

No external library dependencies.

###Visual Studio 2010
- project files included
- requires external [vsyasm.exe](http://yasm.tortall.net/Download.html) in %PATH%
  - run `rundll32 sysdm.cpl,EditEnvironmentVariables` and add PATH to user variables

###Linux
- Simple Makefile included in src/
- Yasm is expected to be in PATH

###Other
- Contact us for support or write an [issue in Github](https://github.com/ultravideo/kvazaar/issues)


##Contributing to Kvazaar

###For version control we try to follow these convetions:

- Master branch always produces a working bitstream (can be decoded with HM).
- Commits for new features and major changes/fixes put to a sensibly named feature branch first and later merged to the master branch.
- Always merge the feature branch to the master branch, not the other way around, with fast-forwarding disabled if necessary. We have found that this differentiates between working and unfinished versions nicely.
- Every commit should at least compile. Producing a working bitstream is nice as well, but not always possible. Features may be temporarily disabled to produce a working bitstream, but remember to re-enbable them before merging to master.


###Testing:

- We do not have a proper testing framework yet. We test mainly by decoding the bitstream with HM and checking that the result matches the encoders own reconstruction.
- You should check that your code encodes at least 600 frames without crashing, that HM can decode the resulting bitstream and produces no hash warnings.
- We would like to have a suite of automatic tests that also check for BD-rate increase and speed decrease in addition to checking that the bitstream is valid. As of yet there is no such suite.
- Compiler should produce no warnings with -W-all. It does now. We are working on fixing that.


###Unit tests:
- There are some unit tests located in the tests directory. We would like to have more.
- The Visual Studio project links the unit tests against the actual .lib file used by the encoder. There is no Makefile as of yet.
- The unit tests use "greatest" unit testing framework. It is included as a submodule, but getting it requires the following commands to be run in the root directory of kvazaar:

        git submodule init
        git submodule update


###Code style:

We try to follow the following conventions:
- ANSI-C/C89.
 - Limited mainly due to poor C support in Visual Studio 2010.
 - // comments allowed and encouraged.
- Follow overall conventions already established in the code.
- Indent by 2 spaces. (no tabs)
- { on the same line for control logic and on the next line for functions
- Reference and deference next to the variable name.
- Variable names in lowered characters with words divided by underscore.
- Maximum line length 79 characters when possible.
- Functions only used inside the module shouldn't be defined in the module header. They can be defined in the beginning of the .c file if necessary.


###Tips for studying HEVC:

- A good first resource is JCTVC-N1002 High Efficiency Video Coding (HEVC) Test Model 12 (HM12) Encoder Description
- Many good articles regarding specific parts of HEVC can be found on IEEE Transactions on Circuits and Systems for Video Technology, Combined issue on High Efficiency Video Coding (HEVC) Standards and Research
- A definitive ansver to a question regarding the bitstream can often be found faster from the HM reference encoder than by reading the specification. 
