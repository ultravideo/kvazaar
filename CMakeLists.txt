cmake_minimum_required(VERSION 3.12)

project(kvazaar
LANGUAGES C CXX
HOMEPAGE_URL https://github.com/ultravideo/kvazaar
DESCRIPTION "An open-source HEVC encoder licensed under 3-clause BSD"
VERSION 2.3.1 )

option(BUILD_SHARED_LIBS "Build using shared kvazaar library" ON)

option(BUILD_TESTS "Build tests" ON)

option(USE_CRYPTO "Use crypto library" OFF)

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
    else()
        set(KVZ_COMPILER_VERSION "VS2022")
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
file(GLOB_RECURSE LIB_SOURCES_STRATEGIES RELATIVE ${PROJECT_SOURCE_DIR} "src/strategies/*.c")

list(APPEND LIB_SOURCES ${LIB_SOURCES_STRATEGIES})

# We also need the libmd5
list(APPEND LIB_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/extras/libmd5.c)

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

if(BUILD_SHARED_LIBS)
  list( APPEND CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_LIBDIR}" "./" "../lib" )
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  add_library(kvazaar SHARED ${LIB_SOURCES})
else()
  add_library(kvazaar STATIC ${LIB_SOURCES})
  if(MSVC) # Fix a linking problem with visual studio when the library is the same name as the binary
    set_target_properties(kvazaar PROPERTIES OUTPUT_NAME libkvazaar)
  endif()
  
endif()

target_include_directories(kvazaar PUBLIC src)
target_include_directories(kvazaar PUBLIC src/extras)
target_include_directories(kvazaar PUBLIC src/strategies)

file(GLOB LIB_SOURCES_STRATEGIES_AVX2 RELATIVE ${PROJECT_SOURCE_DIR} "src/strategies/avx2/*.c")
file(GLOB LIB_SOURCES_STRATEGIES_SSE41 RELATIVE ${PROJECT_SOURCE_DIR} "src/strategies/sse41/*.c")

set(CLI_SOURCES "src/encmain.c" "src/cli.c" "src/cli.h" "src/yuv_io.c" "src/yuv_io.h")

# Add the getopt and pthread for visual studio
if(MSVC) 
  list(APPEND CLI_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/extras/getopt.c ${CMAKE_CURRENT_SOURCE_DIR}/src/threadwrapper/src/pthread.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/threadwrapper/src/semaphore.cpp)
endif()

add_executable(kvazaar-bin ${CLI_SOURCES})

set_target_properties(kvazaar-bin PROPERTIES OUTPUT_NAME kvazaar)
set_target_properties(kvazaar-bin PROPERTIES RUNTIME_OUTPUT_NAME kvazaar)

target_link_libraries(kvazaar-bin PUBLIC kvazaar)

if(MSVC)
  target_include_directories(kvazaar PUBLIC src/threadwrapper/include)
  set_property( SOURCE ${LIB_SOURCES_STRATEGIES_AVX2} APPEND PROPERTY COMPILE_FLAGS "/arch:AVX2" )  
else()
  set_target_properties(kvazaar-bin PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src)
  set_target_properties(kvazaar PROPERTIES SOVERSION "7" VERSION "7.4.0")
  list(APPEND ALLOW_AVX2 "x86_64" "AMD64")
  if(${CMAKE_SYSTEM_PROCESSOR} IN_LIST ALLOW_AVX2) 
    set_property( SOURCE ${LIB_SOURCES_STRATEGIES_AVX2} APPEND PROPERTY COMPILE_FLAGS "-mavx2 -mbmi -mpopcnt -mlzcnt -mbmi2" )
    set_property( SOURCE ${LIB_SOURCES_STRATEGIES_SSE41} APPEND PROPERTY COMPILE_FLAGS "-msse4.1" )    
  endif()
  
  # CYGWIN, MSYS, and MINGW seem to be needing this but in some cases 
  # it might be that the toolset is not properly set, so also use this
  # in cases where we are not sure that it is not needed
  if((NOT MSVC AND NOT LINUX AND NOT APPLE) OR (CYGWIN OR MSYS OR MINGW))
    set(CMAKE_C_FLAGS "-Wa,-muse-unaligned-vector-move ${CMAKE_C_FLAGS}")
  endif()

  
  set(THREADS_PREFER_PTHREAD_FLAG ON)
  find_package(Threads REQUIRED)
  target_link_libraries(kvazaar PUBLIC Threads::Threads)

  include(CheckLibraryExists)

  CHECK_LIBRARY_EXISTS(m sin "" HAVE_LIB_M)

  if (HAVE_LIB_M)
      set(EXTRA_LIBS ${EXTRA_LIBS} m)
  endif (HAVE_LIB_M)
  if (USE_CRYPTO)
      set(EXTRA_LIBS ${EXTRA_LIBS} ${CMAKE_BINARY_DIR}/lib/libcryptopp.a)
  endif ()

  target_link_libraries(kvazaar PUBLIC ${EXTRA_LIBS})
  target_link_libraries(kvazaar-bin PUBLIC ${EXTRA_LIBS} )
endif()


# Source grouping

# Some basic structuring of the files based on previous visual studio project files
file(GLOB SOURCE_GROUP_BITSTREAM RELATIVE ${PROJECT_SOURCE_DIR} "src/encode_coding_tree.*" "src/encoder_state-bitstream.*" "src/nal.*")
file(GLOB SOURCE_GROUP_CABAC RELATIVE ${PROJECT_SOURCE_DIR} "src/bitstream.*" "src/cabac.*" "src/context.*")
file(GLOB SOURCE_GROUP_COMPRESSION RELATIVE ${PROJECT_SOURCE_DIR} "src/search*" "src/rdo.*" "src/fast_coeff*")
file(GLOB SOURCE_GROUP_CONSTRAINT RELATIVE ${PROJECT_SOURCE_DIR} "src/constraint.*" "src/ml_*")
file(GLOB SOURCE_GROUP_CONTROL RELATIVE ${PROJECT_SOURCE_DIR} "src/cfg.*" "src/encoder.*" "src/encoder_state-c*" "src/encoder_state-g*" "src/encoderstate*" "src/gop.*" "src/input_frame_buffer.*" "src/kvazaar*" "src/rate_control.*" "src/mip_data.h")
file(GLOB SOURCE_GROUP_DATA_STRUCTURES RELATIVE ${PROJECT_SOURCE_DIR} "src/cu.*" "src/image.*" "src/imagelist.*" "src/videoframe.*")
file(GLOB SOURCE_GROUP_EXTRAS RELATIVE ${PROJECT_SOURCE_DIR} "src/extras/*.h" "src/extras/*.c")
file(GLOB_RECURSE SOURCE_GROUP_STRATEGIES RELATIVE ${PROJECT_SOURCE_DIR} "src/strategies/*.h" "src/strategies/*.c")
file(GLOB SOURCE_GROUP_RECON RELATIVE ${PROJECT_SOURCE_DIR} "src/alf.*" "src/filter.*" "src/inter.*" "src/intra.*" "src/reshape.*" "src/sao.*" "src/scalinglist.*" "src/tables.*" "src/transform.*")
file(GLOB SOURCE_GROUP_THREADING RELATIVE ${PROJECT_SOURCE_DIR} "src/threadqueue.*" "src/threads.*")
file(GLOB_RECURSE SOURCE_GROUP_THREADWRAPPER RELATIVE ${PROJECT_SOURCE_DIR} "src/threadwrapper/*.cpp" "src/threadwrapper/*.h")
file(GLOB SOURCE_GROUP_TOPLEVEL RELATIVE ${PROJECT_SOURCE_DIR} "src/debug.*" "src/global.h" "src/version.h" "src/kvz_math.h" "src/checkpoint.*")

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
install(TARGETS kvazaar-bin DESTINATION ${CMAKE_INSTALL_BINDIR})
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

  if(NOT "test_external_symbols.sh" IN_LIST XFAIL)
    add_test( NAME test_external_symbols COMMAND ${PROJECT_SOURCE_DIR}/tests/test_external_symbols.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()

  if(NOT "util.sh" IN_LIST XFAIL)
    add_test( NAME util COMMAND ${PROJECT_SOURCE_DIR}/tests/util.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
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

  if(NOT "test_rate_control.sh" IN_LIST XFAIL)
    add_test( NAME test_rate_control COMMAND ${PROJECT_SOURCE_DIR}/tests/test_rate_control.sh WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  endif()
endif()
