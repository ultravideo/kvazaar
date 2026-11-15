cmake_minimum_required(VERSION 3.25)

set(KVAZAAR_VERSION_MAJOR 2)
set(KVAZAAR_VERSION_MINOR 3)
set(KVAZAAR_VERSION_REVISION 2)

project(kvazaar
LANGUAGES C CXX
HOMEPAGE_URL https://github.com/ultravideo/kvazaar
DESCRIPTION "An open-source HEVC encoder licensed under 3-clause BSD"
VERSION ${KVAZAAR_VERSION_MAJOR}.${KVAZAAR_VERSION_MINOR}.${KVAZAAR_VERSION_REVISION})

option(BUILD_SHARED_LIBS "Build using shared kvazaar library" ON)

# Do not fuse floating point multiplications and additions (FMA). The
# encoder decisions (RD costs, rate control, lookahead, temporal filter)
# use floating point, and fusing changes their rounding, so AArch64 (where
# GCC fuses by default) and x86-64 with FMA would code different bitstreams
# than x86-64 without it. MSVC does not fuse under /fp:precise, see also
# global.h.
if(NOT MSVC)
  add_compile_options(-ffp-contract=off)
endif()

option(BUILD_TESTS "Build tests" ON)

option(USE_CRYPTO "Use crypto library" OFF)

# The library contains an encoder for each of these bit depths, see the
# "Bit depths" section below.
set(KVZ_BIT_DEPTHS "8;10;12" CACHE STRING "Bit depths to compile into kvazaar (8, 10 and/or 12)")

option(BUILD_KVAZAAR_BINARY "Build kvazaar-bin" ON) # To build only the lib, useful for iOS-builds

include(GNUInstallDirs) #Helps to define correct distro specific install directories
set(DEFERRED "@")

set(KVAZAAR_INSTALL_LIBDIR "${CMAKE_INSTALL_LIBDIR}" CACHE PATH "kvazaar library install path")
set(KVAZAAR_INSTALL_BINDIR "${CMAKE_INSTALL_BINDIR}" CACHE PATH "kvazaar binary install path")
set(KVAZAAR_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}" CACHE PATH "kvazaar include install path")
set(KVAZAAR_INSTALL_MANDIR "${CMAKE_INSTALL_MANDIR}/man1" CACHE PATH "kvazaar manual page file install path")

# https://www.kitware.com/cmake-and-the-default-build-type/
# Set a default build type if none was specified
set(KVZ_DEFAULT_BUILD_TYPE "RelWithDebInfo")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type specified, setting to '${KVZ_DEFAULT_BUILD_TYPE}'.")
  set(CMAKE_BUILD_TYPE "${KVZ_DEFAULT_BUILD_TYPE}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()


find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    # Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(WARNING "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
    # Check git hash and fetch tag
    execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    RESULT_VARIABLE GIT_HEAD_OK
                    OUTPUT_VARIABLE GIT_HEAD)
    if(GIT_HEAD_OK EQUAL "0")
        string(SUBSTRING ${GIT_HEAD} 0 30 GIT_TAG_LONG)    
        execute_process(COMMAND ${GIT_EXECUTABLE} name-rev --tags --name-only ${GIT_TAG_LONG}
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_TAG_OK
                        OUTPUT_VARIABLE GIT_TAG)
        string(SUBSTRING ${GIT_TAG} 0 9 GIT_TAG_STRIP)
      
        # If tag is not defined, add part of the commit hash to the version
        if(GIT_TAG_OK EQUAL "0" AND GIT_TAG_STRIP STREQUAL "undefined")
          string(SUBSTRING ${GIT_HEAD} 0 7 GIT_TAG_SHORT)
          set(PROJECT_VERSION ${PROJECT_VERSION}-${GIT_TAG_SHORT})
          message(INFO " No tag detected, version changed to ${PROJECT_VERSION}")
        endif()
    endif()    
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/greatest/greatest.h")
    message(WARNING "The submodule greatest was not loaded, some tests may fail")
endif()

# Grab <year>-<month>-<day> timestamp for debug purposes
string(TIMESTAMP CMAKE_BUILD_DATE %Y-%m-%d)

set(KVZ_COMPILER_VERSION "${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}")

if(MSVC)
    if(MSVC_VERSION LESS 1800)
        set(KVZ_COMPILER_VERSION "VS")
    elseif(MSVC_VERSION LESS 1900)
        set(KVZ_COMPILER_VERSION "VS2013")
    elseif(MSVC_VERSION LESS 1910)
        set(KVZ_COMPILER_VERSION "VS2015")
    elseif(MSVC_VERSION LESS 1920)
        set(KVZ_COMPILER_VERSION "VS2017")
    elseif(MSVC_VERSION LESS 1930)
        set(KVZ_COMPILER_VERSION "VS2019")
    elseif(MSVC_VERSION LESS 1950)
        set(KVZ_COMPILER_VERSION "VS2022")
    else()
        set(KVZ_COMPILER_VERSION "VS2026")
    endif()
endif()

# Set compiler info to print at runtime
set(KVZ_COMPILER_STRING "${KVZ_COMPILER_VERSION}")

add_definitions(-DCMAKE_BUILD)

# Set correct pkgconfig libdir variable
if(IS_ABSOLUTE "${CMAKE_INSTALL_LIBDIR}")
  set(KVAZAAR_PC_LIBDIR "${CMAKE_INSTALL_LIBDIR}")
else()
  set(KVAZAAR_PC_LIBDIR "\${exec_prefix}/${CMAKE_INSTALL_LIBDIR}")
endif()

# Set correct pkgconfig include variable
if(IS_ABSOLUTE "${CMAKE_INSTALL_INCLUDEDIR}")
  set(KVAZAAR_PC_INCDIR "${CMAKE_INSTALL_INCLUDEDIR}")
else()
  set(KVAZAAR_PC_INCDIR "\${prefix}/${CMAKE_INSTALL_INCLUDEDIR}")
endif()

# Apply dynamic info to the config files
configure_file("${PROJECT_SOURCE_DIR}/src/kvazaarCMake.pc.in" "${PROJECT_SOURCE_DIR}/src/kvazaar.pc.temp" @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/src/version.h.in" "${PROJECT_SOURCE_DIR}/src/version.h" @ONLY)

# Add all sources in src/ base
file(GLOB LIB_SOURCES RELATIVE ${PROJECT_SOURCE_DIR} "src/*.h" "src/*.c")

# We don't want CLI main in the library
list(REMOVE_ITEM LIB_SOURCES "src/encmain.c" "src/cli.c" "src/cli.h" "src/yuv_io.c" "src/yuv_io.h")

# Add also all the strategies
file(GLOB_RECURSE LIB_SOURCES_STRATEGIES RELATIVE ${PROJECT_SOURCE_DIR} "src/strategies/*.c" "src/strategies/*.cpp")

# SIMD strategies written in C++17 use templates for compile-time specialization
# only (see src/strategies/strategies-cpp.h): no exceptions, RTTI or C++ runtime.
file(GLOB_RECURSE LIB_SOURCES_STRATEGIES_CPP RELATIVE ${PROJECT_SOURCE_DIR} "src/strategies/*.cpp")
if(MSVC)
  set_property(SOURCE ${LIB_SOURCES_STRATEGIES_CPP} APPEND PROPERTY COMPILE_OPTIONS "/std:c++17")
else()
  set_property(SOURCE ${LIB_SOURCES_STRATEGIES_CPP} APPEND PROPERTY COMPILE_OPTIONS "-std=c++17;-fno-exceptions;-fno-rtti")
endif()

list(APPEND LIB_SOURCES ${LIB_SOURCES_STRATEGIES})

# The AArch64 NEON strategies (src/strategies/neon) are picked up by the glob
# above and need no extra compiler flags: NEON is part of the AArch64 baseline.

# We also need the libmd5
list(APPEND LIB_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/extras/libmd5.c)

# Sources that are compiled for every bit depth. The platform support code
# and the C++ encryption code added below do not depend on the bit depth
# and are only compiled once.
set(LIB_SOURCES_BIT_DEPTH ${LIB_SOURCES})

add_definitions(-DKVZ_DLL_EXPORTS)

if(BUILD_SHARED_LIBS)
  add_definitions(-DPIC)
endif()

# For visual studio / windows we also need our own pthread implementation and getopt
if(MSVC)
  list(APPEND LIB_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/extras/getopt.c ${CMAKE_CURRENT_SOURCE_DIR}/src/threadwrapper/src/pthread.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/threadwrapper/src/semaphore.cpp)
  add_definitions(-DWIN32_LEAN_AND_MEAN -D_WIN32 -DWIN32 -DWIN64)
endif()


if (USE_CRYPTO)
    if(BUILD_TESTS)
        message(WARNING "Crypto++ is not compatible with the tests, disabling tests")
        set(BUILD_TESTS OFF)
    endif()
    if (BUILD_SHARED_LIBS)
        message(WARNING "Crypto++ is not compatible with shared libraries, disabling shared libraries")
        set(BUILD_SHARED_LIBS OFF)
    endif ()
    include(FetchContent)

    add_definitions(-DKVZ_SEL_ENCRYPTION)
    list(APPEND LIB_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/extras/crypto.cpp)
    list(APPEND CLI_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/extras/crypto.cpp)

    if (NOT CRYPTOPP_FOUND)
        message(STATUS "Fetching and building Crypto++")

        # CMake for Crypto++
        FetchContent_Declare(
                cryptopp-cmake
                GIT_REPOSITORY https://github.com/abdes/cryptopp-cmake.git
                GIT_TAG        43367a9cef6576b34179427a31a619802205406e
        )

        set(CRYPTOPP_INSTALL OFF CACHE BOOL "" FORCE) # we don't want to install Crypto++
        set(CRYPTOPP_BUILD_TESTING OFF CACHE BOOL "" FORCE)
        set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)

        FetchContent_MakeAvailable(cryptopp-cmake)

        unset(BUILD_SHARED_LIBS) # unset so it does not affect other projects

        # copy lib binary so it is found later
        file(GLOB CRYPTOPP_BIN "${cryptopp-cmake_BINARY_DIR}/cryptopp/cryptopp.*")
        file(COPY ${CRYPTOPP_BIN} DESTINATION ${CMAKE_BINARY_DIR}/lib/)
        file(GLOB CRYPTOPP_BIN "${cryptopp-cmake_BINARY_DIR}/cryptopp/libcryptopp.*")
        file(COPY ${CRYPTOPP_BIN} DESTINATION ${CMAKE_BINARY_DIR}/lib/)
    endif ()


endif ()

# The library sources are compiled into an object library, which the tests
# also use (see tests/CMakeLists.txt). Like with autotools, only the public
# API (KVZ_PUBLIC) is visible outside of the shared library, so the compiler
# can inline the internal functions and call them directly instead of through
# the PLT.
add_library(kvazaar_objects OBJECT ${LIB_SOURCES})
set_target_properties(kvazaar_objects PROPERTIES
  POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS}
  C_VISIBILITY_PRESET hidden
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON)

if(BUILD_SHARED_LIBS)
  list( APPEND CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_LIBDIR}" "./" "../lib" )
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  add_library(kvazaar SHARED)
else()
  add_library(kvazaar STATIC)
  if(MSVC) # Fix a linking problem with visual studio when the library is the same name as the binary
    set_target_properties(kvazaar PROPERTIES OUTPUT_NAME libkvazaar)
  endif()
  
endif()

# The C++ strategy files need no C++ runtime, so link the library as C to avoid
# a libstdc++ dependency. MSVC and the Crypto++ build do use the runtime.
# The library has no sources of its own, so set the language in both cases.
# Linking as C also uses the C flags, such as -fsanitize=address.
if(NOT MSVC AND NOT USE_CRYPTO)
  set_target_properties(kvazaar PROPERTIES LINKER_LANGUAGE C)
else()
  set_target_properties(kvazaar PROPERTIES LINKER_LANGUAGE CXX)
endif()

# Bit depths
#
# The lowest bit depth of KVZ_BIT_DEPTHS is the default one, compiled as
# usual. The others are compiled into object libraries whose external
# symbols are renamed by force including src/multilib_rename.h (e.g.
# kvz_api_get becomes kvz_api_get_10bit), and linked into the same library.
# kvz_api_get of the default bit depth returns the API of the others.
string(REPLACE "," ";" KVZ_BIT_DEPTHS "${KVZ_BIT_DEPTHS}")
if(CMAKE_C_FLAGS MATCHES "-DKVZ_BIT_DEPTH=([0-9]+)")
  # Old way of selecting a single bit depth.
  set(KVZ_BIT_DEPTHS ${CMAKE_MATCH_1})
  message(WARNING "-DKVZ_BIT_DEPTH in CMAKE_C_FLAGS is deprecated, use -DKVZ_BIT_DEPTHS=${KVZ_BIT_DEPTHS} instead.")
endif()
# The bit depth is given to the C and the C++ sources by the targets below.
string(REGEX REPLACE "-DKVZ_BIT_DEPTH=[0-9]+" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
string(REGEX REPLACE "-DKVZ_BIT_DEPTH=[0-9]+" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
list(REMOVE_DUPLICATES KVZ_BIT_DEPTHS)
list(SORT KVZ_BIT_DEPTHS COMPARE NATURAL)
foreach(depth IN LISTS KVZ_BIT_DEPTHS)
  if(NOT depth MATCHES "^(8|10|12)$")
    message(FATAL_ERROR "Unsupported bit depth in KVZ_BIT_DEPTHS: ${depth}. Supported bit depths are 8, 10 and 12.")
  endif()
endforeach()
list(LENGTH KVZ_BIT_DEPTHS num_bit_depths)
if(num_bit_depths EQUAL 0)
  message(FATAL_ERROR "KVZ_BIT_DEPTHS is empty.")
endif()
list(GET KVZ_BIT_DEPTHS 0 KVZ_DEFAULT_BIT_DEPTH)
set(KVZ_EXTRA_BIT_DEPTHS ${KVZ_BIT_DEPTHS})
list(REMOVE_AT KVZ_EXTRA_BIT_DEPTHS 0)
message(STATUS "Bit depths: ${KVZ_BIT_DEPTHS}")

# Only the objects linked directly are added to a target, so kvazaar gets
# the objects and the targets linking to kvazaar the usage requirements.
target_link_libraries(kvazaar PUBLIC kvazaar_objects)

target_compile_definitions(kvazaar_objects PUBLIC KVZ_BIT_DEPTH=${KVZ_DEFAULT_BIT_DEPTH})
target_compile_definitions(kvazaar_objects PUBLIC KVZ_VERSION_MAJOR=${KVAZAAR_VERSION_MAJOR})
target_compile_definitions(kvazaar_objects PUBLIC KVZ_VERSION_MINOR=${KVAZAAR_VERSION_MINOR})
target_compile_definitions(kvazaar_objects PUBLIC KVZ_VERSION_REVISION=${KVAZAAR_VERSION_REVISION})

target_include_directories(kvazaar_objects PUBLIC src)
target_include_directories(kvazaar_objects PUBLIC src/extras)
target_include_directories(kvazaar_objects PUBLIC src/strategies)

file(GLOB LIB_SOURCES_STRATEGIES_AVX2 RELATIVE ${PROJECT_SOURCE_DIR} "src/strategies/avx2/*.c" "src/strategies/avx2/*.cpp")
file(GLOB LIB_SOURCES_STRATEGIES_SSE41 RELATIVE ${PROJECT_SOURCE_DIR} "src/strategies/sse41/*.c")

set(CLI_SOURCES "src/encmain.c" "src/cli.c" "src/cli.h" "src/yuv_io.c" "src/yuv_io.h")

# Add the getopt and pthread for visual studio
if(MSVC) 
  list(APPEND CLI_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/extras/getopt.c ${CMAKE_CURRENT_SOURCE_DIR}/src/threadwrapper/src/pthread.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/threadwrapper/src/semaphore.cpp)
endif()

if(BUILD_KVAZAAR_BINARY)
  add_executable(kvazaar-bin ${CLI_SOURCES})
  if(NOT MSVC AND NOT USE_CRYPTO)
    set_target_properties(kvazaar-bin PROPERTIES LINKER_LANGUAGE C)
  endif()

  if(MSVC AND BUILD_SHARED_LIBS)
    set_target_properties(kvazaar PROPERTIES
      OUTPUT_NAME kvazaar
      PDB_NAME kvazaar-lib
    )

    set_target_properties(kvazaar-bin PROPERTIES
      OUTPUT_NAME kvazaar
      PDB_NAME kvazaar-bin
    )
  else()
    set_target_properties(kvazaar-bin PROPERTIES
      OUTPUT_NAME kvazaar
    )
  endif()
  set_target_properties(kvazaar-bin PROPERTIES RUNTIME_OUTPUT_NAME kvazaar)

  target_link_libraries(kvazaar-bin PUBLIC kvazaar)
endif()

list(APPEND ALLOW_AVX2 "x86_64" "AMD64")

if(MSVC)
  target_include_directories(kvazaar_objects PUBLIC src/threadwrapper/include)
  if(${CMAKE_SYSTEM_PROCESSOR} IN_LIST ALLOW_AVX2) 
    set_property( SOURCE ${LIB_SOURCES_STRATEGIES_AVX2} APPEND PROPERTY COMPILE_FLAGS "/arch:AVX2" )  
  endif()
else()
  if(BUILD_KVAZAAR_BINARY)
    set_target_properties(kvazaar-bin PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src)
  endif()
  set_target_properties(kvazaar PROPERTIES SOVERSION "7" VERSION "7.6.0")
  if(${CMAKE_SYSTEM_PROCESSOR} IN_LIST ALLOW_AVX2) 
    set_property( SOURCE ${LIB_SOURCES_STRATEGIES_AVX2} APPEND PROPERTY COMPILE_FLAGS "-mavx2 -mbmi -mpopcnt -mlzcnt -mbmi2" )
    set_property( SOURCE ${LIB_SOURCES_STRATEGIES_SSE41} APPEND PROPERTY COMPILE_FLAGS "-msse4.1" )    
  endif()
  
  # CYGWIN, MSYS, and MINGW seem to be needing this but in some cases 
  # it might be that the toolset is not properly set, so also use this
  # in cases where we are not sure that it is not needed
  if((NOT MSVC AND NOT LINUX AND NOT APPLE AND NOT ANDROID AND NOT EMSCRIPTEN AND NOT BSD) OR (CYGWIN OR MSYS OR MINGW))
    set(CMAKE_C_FLAGS "-Wa,-muse-unaligned-vector-move ${CMAKE_C_FLAGS}")
  endif()

  if(NOT ANDROID AND NOT EMSCRIPTEN)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
    target_link_libraries(kvazaar_objects PUBLIC Threads::Threads)
  endif()

  include(CheckLibraryExists)

  CHECK_LIBRARY_EXISTS(m sin "" HAVE_LIB_M)

  if (HAVE_LIB_M)
      set(EXTRA_LIBS ${EXTRA_LIBS} m)
  endif (HAVE_LIB_M)
  if (USE_CRYPTO)
      set(EXTRA_LIBS ${EXTRA_LIBS} ${CMAKE_BINARY_DIR}/lib/libcryptopp.a)
  endif ()

  target_link_libraries(kvazaar_objects PUBLIC ${EXTRA_LIBS})
  if(BUILD_KVAZAAR_BINARY)
    target_link_libraries(kvazaar-bin PUBLIC ${EXTRA_LIBS} )
  endif()
endif()


foreach(depth IN LISTS KVZ_EXTRA_BIT_DEPTHS)
  set(depth_target kvazaar_${depth}bit)
  add_library(${depth_target} OBJECT ${LIB_SOURCES_BIT_DEPTH})
  set_target_properties(${depth_target} PROPERTIES
    POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS}
    C_VISIBILITY_PRESET hidden
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON)
  target_compile_definitions(${depth_target} PRIVATE
    KVZ_BIT_DEPTH=${depth}
    KVZ_VERSION_MAJOR=${KVAZAAR_VERSION_MAJOR}
    KVZ_VERSION_MINOR=${KVAZAAR_VERSION_MINOR}
    KVZ_VERSION_REVISION=${KVAZAAR_VERSION_REVISION})
  target_include_directories(${depth_target} PRIVATE src src/extras src/strategies)
  if(MSVC)
    target_include_directories(${depth_target} PRIVATE src/threadwrapper/include)
    target_compile_options(${depth_target} PRIVATE "/FI${PROJECT_SOURCE_DIR}/src/multilib_rename.h")
  else()
    target_compile_options(${depth_target} PRIVATE -include "${PROJECT_SOURCE_DIR}/src/multilib_rename.h")
    if(TARGET Threads::Threads)
      target_link_libraries(${depth_target} PRIVATE Threads::Threads)
    endif()
  endif()

  target_compile_definitions(kvazaar_objects PRIVATE KVZ_MULTILIB_${depth}BIT)
  target_sources(kvazaar PRIVATE $<TARGET_OBJECTS:${depth_target}>)
  list(APPEND KVZ_EXTRA_BIT_DEPTH_TARGETS ${depth_target})
endforeach()


# Source grouping

# Some basic structuring of the files based on previous visual studio project files
file(GLOB SOURCE_GROUP_BITSTREAM RELATIVE ${PROJECT_SOURCE_DIR} "src/encode_coding_tree.*" "src/encoder_state-bitstream.*" "src/nal.*")
file(GLOB SOURCE_GROUP_CABAC RELATIVE ${PROJECT_SOURCE_DIR} "src/bitstream.*" "src/cabac.*" "src/context.*")
file(GLOB SOURCE_GROUP_COMPRESSION RELATIVE ${PROJECT_SOURCE_DIR} "src/search*" "src/rdo.*" "src/fast_coeff*" "src/work_tree_pool.*")
file(GLOB SOURCE_GROUP_CONSTRAINT RELATIVE ${PROJECT_SOURCE_DIR} "src/constraint.*" "src/ml_*")
file(GLOB SOURCE_GROUP_CONTROL RELATIVE ${PROJECT_SOURCE_DIR} "src/cfg.*" "src/encoder.*" "src/encoder_state-c*" "src/encoder_state-g*" "src/encoderstate*" "src/gop.*" "src/input_frame_buffer.*" "src/lookahead.*" "src/temporal_filter.*" "src/kvazaar*" "src/rate_control.*" "src/mip_data.h")
file(GLOB SOURCE_GROUP_DATA_STRUCTURES RELATIVE ${PROJECT_SOURCE_DIR} "src/cu.*" "src/image.*" "src/imagelist.*" "src/videoframe.*")
file(GLOB SOURCE_GROUP_EXTRAS RELATIVE ${PROJECT_SOURCE_DIR} "src/extras/*.h" "src/extras/*.c")
file(GLOB_RECURSE SOURCE_GROUP_STRATEGIES RELATIVE ${PROJECT_SOURCE_DIR} "src/strategies/*.h" "src/strategies/*.c" "src/strategies/*.cpp")
file(GLOB SOURCE_GROUP_RECON RELATIVE ${PROJECT_SOURCE_DIR} "src/alf.*" "src/filter.*" "src/inter.*" "src/intra.*" "src/reshape.*" "src/sao.*" "src/scalinglist.*" "src/tables.*" "src/transform.*")
file(GLOB SOURCE_GROUP_THREADING RELATIVE ${PROJECT_SOURCE_DIR} "src/threadqueue.*" "src/threads.*")
file(GLOB_RECURSE SOURCE_GROUP_THREADWRAPPER RELATIVE ${PROJECT_SOURCE_DIR} "src/threadwrapper/*.cpp" "src/threadwrapper/*.h")
file(GLOB SOURCE_GROUP_TOPLEVEL RELATIVE ${PROJECT_SOURCE_DIR} "src/debug.*" "src/global.h" "src/version.h" "src/kvz_math.h" "src/checkpoint.*" "src/multilib_rename.h")

source_group( "Bitstream" FILES  ${SOURCE_GROUP_BITSTREAM})
source_group( "CABAC" FILES ${SOURCE_GROUP_CABAC})
source_group( "Compression" FILES ${SOURCE_GROUP_COMPRESSION})
source_group( "Constraint" FILES ${SOURCE_GROUP_CONSTRAINT})
source_group( "Control" FILES ${SOURCE_GROUP_CONTROL})
source_group( "Data Structures" FILES ${SOURCE_GROUP_DATA_STRUCTURES})
source_group( "Extras" FILES ${SOURCE_GROUP_EXTRAS})

# Handle the strategies directory structure better in visual studio
if(MSVC)
  foreach(source IN LISTS SOURCE_GROUP_STRATEGIES)
      get_filename_component(source_path "${source}" PATH)
      string(REPLACE "src/" "" source_path_msvc "${source_path}")
      string(REPLACE "/" "\\" source_path_msvc "${source_path_msvc}")
      source_group("Optimization\\${source_path_msvc}" FILES "${source}")
  endforeach()
else()
  source_group( "Optimization" FILES ${SOURCE_GROUP_STRATEGIES})
endif()
source_group( "Optimization" FILES "src/strategyselector.c" "src/strategyselector.h")

source_group( "Reconstruction" FILES ${SOURCE_GROUP_RECON})
source_group( "Threading" FILES ${SOURCE_GROUP_THREADING})
source_group( "Threadwrapper" FILES ${SOURCE_GROUP_THREADWRAPPER})
source_group( "" FILES ${SOURCE_GROUP_TOPLEVEL})

# INSTALL

# ToDo: make configurable
install(CODE "configure_file(\"${PROJECT_SOURCE_DIR}/src/kvazaar.pc.temp\" \"${PROJECT_SOURCE_DIR}/src/kvazaar.pc\" @ONLY)")

install(FILES ${PROJECT_SOURCE_DIR}/src/kvazaar.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
if(BUILD_KVAZAAR_BINARY)
  install(TARGETS kvazaar-bin DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
install(TARGETS kvazaar
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${PROJECT_SOURCE_DIR}/src/kvazaar.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${PROJECT_SOURCE_DIR}/doc/kvazaar.1 DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)

IF(UNIX)
# DIST

set(GIT_LS_TREE_OK "1")

# By default grab the list of files in the git repo
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    execute_process(COMMAND ${GIT_EXECUTABLE} ls-tree --name-only -r HEAD
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    RESULT_VARIABLE GIT_LS_TREE_OK
                    OUTPUT_VARIABLE GIT_LS_TREE)
    if(GIT_LS_TREE_OK EQUAL "0")
        string(REGEX REPLACE "\n" ";" GIT_LS_TREE "${GIT_LS_TREE}")
        string(REGEX REPLACE "\r" "" GIT_LS_TREE "${GIT_LS_TREE}")
        list(APPEND DIST_SOURCES ${GIT_LS_TREE})
    endif()
endif()
if(NOT GIT_LS_TREE_OK EQUAL "0")
    file(GLOB_RECURSE DIST_SOURCES RELATIVE ${PROJECT_SOURCE_DIR} "src/*.c" "src/*.h" "tests/*.sh" "tools/*.sh" "tools/*.py" ".github/*.yml" "src/*.in" "placeholder.txt" "CMakeLists.txt" "doc/*" "examples/*" "rdcost-weight-tool/*" "greatest/*.h" "greatest/*.md")
    list(APPEND DIST_SOURCES ".clang-format" ".gitignore" ".gitmodules" "tests/tsan_suppressions.txt" ".travis-install.bash" "CREDITS" "Dockerfile" "docs.doxy" ".gitlab-ci.yml" "LICENSE" "LICENSE.EXT.greatest" "README.md")
endif()

add_custom_target(dist
                COMMAND echo \"Writing log to ${PROJECT_SOURCE_DIR}/dist.log\"
                 && tar -zcvf "${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz" --transform 's,^,${PROJECT_NAME}-${PROJECT_VERSION}/,' -- ${DIST_SOURCES} > dist.log 2>&1 || { echo \"\\033[0;31mfailed to pack ${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz, check ${PROJECT_SOURCE_DIR}/dist.log.\\033\[m\"$<SEMICOLON> exit 1$<SEMICOLON> }
                COMMENT "Make distribution ${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz"
                WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
                BYPRODUCTS ${CMAKE_SOURCE_DIR}/${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz
            )

# DISTCHECK

set(TEMP_DISTCHECK_DIR "_distcheck")

add_custom_target(distcheck
                  COMMAND echo \"Writing log to ${PROJECT_SOURCE_DIR}/distcheck.log\"
                  && cd ${PROJECT_SOURCE_DIR}
                  && mkdir -p ${TEMP_DISTCHECK_DIR}
                  && cd ${TEMP_DISTCHECK_DIR}
                  && tar -zxf ${CMAKE_SOURCE_DIR}/${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz > ${PROJECT_SOURCE_DIR}/distcheck.log || { echo \"\\033[0;31mfailed to unpack ${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz.\\033\[m\"$<SEMICOLON> exit 1$<SEMICOLON> }
                  && echo \"\\033[0;32mFile unpack ok\\033[m\"
                  && cd ${PROJECT_NAME}-${PROJECT_VERSION}
                  && mkdir -p build
                  && cd build                  
                  && cmake -DCMAKE_INSTALL_PREFIX=./ -DBUILD_SHARED_LIBS=OFF -G "Unix Makefiles" .. >> ${PROJECT_SOURCE_DIR}/distcheck.log || { echo \"\\033[0;31mcmake failed to configure.\\033[m\"$<SEMICOLON> exit 1$<SEMICOLON> }
                  && echo \"\\033[0;32mCMake configure ok\\033[m\"
                  && make -j >> ${PROJECT_SOURCE_DIR}/distcheck.log || { echo \"\\033[0;31mmake failed.\\033[m\"$<SEMICOLON> exit 1$<SEMICOLON> }
                  && echo \"\\033[0;32mMake ok\\033[m\"
                  # Full tests might be too demanding to run, enable with parameter?
                  #&& make test || (echo \"\\e[0;31mmake test failed.\\033[m\" && false)
                  && tests/kvazaar_tests >> ${PROJECT_SOURCE_DIR}/distcheck.log 2>&1 || { echo \"\\033[0;31mtests failed.\\033[m\"$<SEMICOLON> exit 1$<SEMICOLON> }
                  && echo \"\\033[0;32mTests ok\\033[m\"
                  && make install >> ${PROJECT_SOURCE_DIR}/distcheck.log || { echo \"\\033[0;31mmake install failed.\\033[m\"$<SEMICOLON> exit 1$<SEMICOLON> }
                  && echo \"\\033[0;32mInstall ok\\033[m\"
                  && bin/kvzaar --help >> ${PROJECT_SOURCE_DIR}/distcheck.log || { echo \"\\033[0;31mkvazaar binary failed to run.\\033[m\"$<SEMICOLON> exit 1$<SEMICOLON> }
                  && echo \"\\033[0;32mbin/kvazaar ok\\033[m\"
                  && make clean >> ${PROJECT_SOURCE_DIR}/distcheck.log || { echo \"\\033[0;31mmake clean failed.\\033[m\"$<SEMICOLON> exit 1$<SEMICOLON> }
                  && echo \"\\033[0;32mmake clean ok\\033[m\"
                  && cd ${PROJECT_SOURCE_DIR}
                  && rm -rf "${PROJECT_SOURCE_DIR}/${TEMP_DISTCHECK_DIR}"
                  && echo \"\\033[0;32m==============================================================\\033[m\"
                  && echo \"\\033[0;32m${PROJECT_NAME}-${PROJECT_VERSION} archives ready for distribution:\\033[m\"
                  && echo \"\\033[0;32m${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz\\033[m\"
                  && echo \"\\033[0;32m==============================================================\\033[m\"
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                  DEPENDS ${CMAKE_SOURCE_DIR}/${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz
                  COMMENT "Checking ${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz.."
                  )
endif() #Unix

# TESTS
enable_testing()

if(MSVC OR MINGW OR MSYS)
  if(BUILD_SHARED_LIBS)
    set(BUILD_TESTS OFF)
    message(INFO " Disable test building, fails in MSVC/MINGW/MSYS2 when building shared binaries")
  endif()
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/greatest/greatest.h" AND BUILD_TESTS)
  add_subdirectory( "tests/" )
  add_test( NAME Test_kvazaar COMMAND kvazaar_tests )
endif()

if(NOT DEFINED MSVC)
  list(APPEND XFAIL "off")
  if(DEFINED ENV{XFAIL_TESTS})
    list(APPEND XFAIL $ENV{XFAIL_TESTS})
  endif()

  if(NOT "test_tools.sh" IN_LIST XFAIL)
    add_test( NAME test_tools COMMAND ${PROJECT_SOURCE_DIR}/tests/test_tools.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()

  if(NOT "test_smp.sh" IN_LIST XFAIL)
    add_test( NAME test_smp COMMAND ${PROJECT_SOURCE_DIR}/tests/test_smp.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()

  if(NOT "test_pu_depth_constraints.sh" IN_LIST XFAIL)
    add_test( NAME test_pu_depth_constraints COMMAND ${PROJECT_SOURCE_DIR}/tests/test_pu_depth_constraints.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()

  if(NOT "test_mv_constraint.sh" IN_LIST XFAIL)
    add_test( NAME test_mv_constraint COMMAND ${PROJECT_SOURCE_DIR}/tests/test_mv_constraint.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()

  if(NOT "test_interlace.sh" IN_LIST XFAIL)
    add_test( NAME test_interlace COMMAND ${PROJECT_SOURCE_DIR}/tests/test_interlace.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()

  if(NOT "test_gop.sh" IN_LIST XFAIL)
    add_test( NAME test_gop COMMAND ${PROJECT_SOURCE_DIR}/tests/test_gop.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()

  if(NOT "test_owf_wpp_tiles.sh" IN_LIST XFAIL)
    add_test( NAME test_owf_wpp_tiles COMMAND ${PROJECT_SOURCE_DIR}/tests/test_owf_wpp_tiles.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()

  if(NOT "test_weird_shapes.sh" IN_LIST XFAIL)
    add_test( NAME test_weird_shapes COMMAND ${PROJECT_SOURCE_DIR}/tests/test_weird_shapes.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()

  # Checks the external symbols of the library objects, so it needs a static
  # library: the shared library only exports the API.
  if(BUILD_SHARED_LIBS AND TARGET kvazaar_tests_lib)
    set(KVZ_STATIC_LIBRARY_TARGET kvazaar_tests_lib)
  elseif(NOT BUILD_SHARED_LIBS)
    set(KVZ_STATIC_LIBRARY_TARGET kvazaar)
  endif()
  if(NOT "test_external_symbols.sh" IN_LIST XFAIL AND KVZ_STATIC_LIBRARY_TARGET)
    add_test( NAME test_external_symbols COMMAND ${PROJECT_SOURCE_DIR}/tests/test_external_symbols.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
    set_tests_properties(test_external_symbols PROPERTIES
      ENVIRONMENT "KVZ_STATIC_LIBRARY=$<TARGET_FILE:${KVZ_STATIC_LIBRARY_TARGET}>")
  endif()

  if(NOT "test_invalid_input.sh" IN_LIST XFAIL)
    add_test( NAME test_invalid_input COMMAND ${PROJECT_SOURCE_DIR}/tests/test_invalid_input.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()

  if(NOT "test_slices.sh" IN_LIST XFAIL)
    add_test( NAME test_slices COMMAND ${PROJECT_SOURCE_DIR}/tests/test_slices.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()

  if(NOT "test_intra.sh" IN_LIST XFAIL)
    add_test( NAME test_intra COMMAND ${PROJECT_SOURCE_DIR}/tests/test_intra.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()

  if(NOT "test_lookahead.sh" IN_LIST XFAIL)
    add_test( NAME test_lookahead COMMAND ${PROJECT_SOURCE_DIR}/tests/test_lookahead.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()

  if(NOT "test_rate_control.sh" IN_LIST XFAIL)
    add_test( NAME test_rate_control COMMAND ${PROJECT_SOURCE_DIR}/tests/test_rate_control.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()
  
  if(NOT "test_chroma_subsampling.sh" IN_LIST XFAIL)
    add_test( NAME test_chroma_subsampling COMMAND ${PROJECT_SOURCE_DIR}/tests/test_chroma_subsampling.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()
  if(NOT "test_chroma_subsampling_422.sh" IN_LIST XFAIL)
    add_test( NAME test_chroma_subsampling_422 COMMAND ${PROJECT_SOURCE_DIR}/tests/test_chroma_subsampling_422.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()
  if(NOT "test_multilib.sh" IN_LIST XFAIL)
    add_test( NAME test_multilib COMMAND ${PROJECT_SOURCE_DIR}/tests/test_multilib.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()
endif()
