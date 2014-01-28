#Kvazaar
An open-source HEVC encoder

TODO: add description

#Using Kvazaar
Currently most of the features are turned on/off from the code on compile time, but they are
meant to be user configurable later.

    kvazaar -i <INPUT_YUV> -w <WIDTH> -h <HEIGHT> -o <OUTPUT.BIN> -n <NUMBER_OF_FRAMES> -q <QP>

eg. `kvazaar -i BQMall_832x480_60.yuv -w 832 -h 480 -o out.bin -n 600 -q 32`

TODO: add examples

##Compiling Kvazaar

If you have trouble regarding compiling the source code, please make an issue about in Github. Others might encounter the same problem and there is probably much to improve in the build process. We want to make this as simple as possible.

- Visual Studio 2010 project included.
- Simple linux Makefile included.
- Yasm used to compile assembly.
- No external library dependencies.

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
- There are some unit tests located in the tests directory. We would like to have more.
- Compiler should produce no warnings with -W-all. It does now. We are working on fixing that.



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
