cmake_minimum_required(VERSION 3.12)

project(kvazaar
LANGUAGES C CXX
HOMEPAGE_URL https://github.com/ultravideo/kvazaar
DESCRIPTION "An open-source VVC encoder licensed under 3-clause BSD"
VERSION 0.2.3 )

# Apply dynamic info to the config files
configure_file("${PROJECT_SOURCE_DIR}/src/kvazaar.pc.in" "${PROJECT_SOURCE_DIR}/src/kvazaar.pc" @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/src/version.h.in" "${PROJECT_SOURCE_DIR}/src/version.h" @ONLY)

option(USE_SHARED_LIB "Build using shared kvazaar library" ON)


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
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/greatest/greatest.h")
    message(WARNING "The submodule greatest was not loaded, some tests may fail")
endif()

# Add all sources in src/ base
file(GLOB LIB_SOURCES RELATIVE ${PROJECT_SOURCE_DIR} "src/*.h" "src/*.c")

# We don't want CLI main in the library
list(REMOVE_ITEM LIB_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/encmain.c" "${CMAKE_CURRENT_SOURCE_DIR}/src/cli.c" "${CMAKE_CURRENT_SOURCE_DIR}/src/yuv_io.c")

# Add also all the strategies
file(GLOB_RECURSE LIB_SOURCES_STRATEGIES RELATIVE ${PROJECT_SOURCE_DIR} "src/strategies/*.c")

# ToDo: do something with encode_coding_tree-avx2, currently not converted to VVC
list(REMOVE_ITEM LIB_SOURCES_STRATEGIES "${CMAKE_CURRENT_SOURCE_DIR}/src/strategies/avx2/encode_coding_tree-avx2.c")

list(APPEND LIB_SOURCES ${LIB_SOURCES_STRATEGIES})

# We also need the libmd5
list(APPEND LIB_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/extras/libmd5.c)

if(NOT USE_SHARED_LIB)
  add_definitions(-DPIC -DKVZ_DLL_EXPORTS)
endif()

# For visual studio / windows we also need our own pthread implementation and getopt
if(MSVC)
  list(APPEND LIB_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/extras/getopt.c ${CMAKE_CURRENT_SOURCE_DIR}/src/threadwrapper/src/pthread.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/threadwrapper/src/semaphore.cpp)
  add_definitions(-DWIN32_LEAN_AND_MEAN -D_WIN32 -DWIN32 -DWIN64)
endif()

# ToDo: allow compiling on other than x86
add_definitions(-DCOMPILE_INTEL)

if(USE_SHARED_LIB)
  add_library(kvazaar SHARED ${LIB_SOURCES})
else()
  add_library(kvazaar STATIC ${LIB_SOURCES})
endif()

target_include_directories(kvazaar PUBLIC src)
target_include_directories(kvazaar PUBLIC src/extras)
target_include_directories(kvazaar PUBLIC src/strategies)

file(GLOB LIB_SOURCES_STRATEGIES_AVX2 RELATIVE ${PROJECT_SOURCE_DIR} "src/strategies/avx2/*.c")

set(CLI_SOURCES "src/encmain.c" "src/cli.c" "src/yuv_io.c")

# Add the getopt and pthread for visual studio
if(MSVC) 
  list(APPEND CLI_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/extras/getopt.c ${CMAKE_CURRENT_SOURCE_DIR}/src/threadwrapper/src/pthread.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/threadwrapper/src/semaphore.cpp)
endif()

add_executable(kvazaar-bin ${CLI_SOURCES})

set_target_properties(kvazaar-bin PROPERTIES OUTPUT_NAME kvazaar)
set_target_properties(kvazaar-bin PROPERTIES RUNTIME_OUTPUT_NAME kvazaar)

target_link_libraries(kvazaar-bin PUBLIC kvazaar)

if(MSVC)
  target_include_directories(kvazaar-bin PUBLIC src/threadwrapper/include)
  set_property( SOURCE ${LIB_SOURCES_STRATEGIES_AVX2} APPEND PROPERTY COMPILE_FLAGS "/arch:AVX2" )  
else()
  set_property( SOURCE ${LIB_SOURCES_STRATEGIES_AVX2} APPEND PROPERTY COMPILE_FLAGS "-mavx2 -mbmi -mpopcnt -mlzcnt -mbmi2" )
  set(THREADS_PREFER_PTHREAD_FLAG ON)
  find_package(Threads REQUIRED)
  target_link_libraries(kvazaar PUBLIC Threads::Threads)

  include(CheckLibraryExists)

  CHECK_LIBRARY_EXISTS(m sin "" HAVE_LIB_M)

  if (HAVE_LIB_M)
      set(EXTRA_LIBS ${EXTRA_LIBS} m)
  endif (HAVE_LIB_M)

  target_link_libraries(kvazaar-bin PUBLIC ${EXTRA_LIBS})
endif()

# INSTALL

# ToDo: make configurable

install(FILES ${PROJECT_SOURCE_DIR}/src/kvazaar.pc DESTINATION ${CMAKE_INSTALL_PREFIX}/share/pkgconfig)
install(TARGETS kvazaar-bin DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
install(TARGETS libkvazaar DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
if(USE_SHARED_LIB) # Just add the lib to the bin directory for now
  install(TARGETS kvazaar DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
endif()
install(FILES ${PROJECT_SOURCE_DIR}/src/kvazaar.h DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
install(FILES ${PROJECT_SOURCE_DIR}/doc/kvazaar.1 DESTINATION ${CMAKE_INSTALL_PREFIX}/share/man)

# DIST

file(GLOB_RECURSE DIST_SOURCES RELATIVE ${PROJECT_SOURCE_DIR} "*.c" "*.h" "*.sh" "*.py" "*.md" "*.yml" "*.in" "placeholder.txt" "CMakeLists.txt" "doc/*" "greatest/*")
list(APPEND DIST_SOURCES ".clang-format" ".gitignore" ".gitmodules" "tests/tsan_suppressions.txt" ".travis-install.bash" "CREDITS" "Dockerfile" "docs.doxy" "LICENSE" "LICENSE.EXT.greatest")

add_custom_target(dist
                COMMAND tar -zcvf "${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz" --transform 's,^,${PROJECT_NAME}-${PROJECT_VERSION}/,' -- ${DIST_SOURCES}
                COMMENT "Make distribution ${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz"
                WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
            )

# TESTS
if(EXISTS "${PROJECT_SOURCE_DIR}/greatest/greatest.h")
  add_subdirectory( "tests/" )
  add_test( NAME Test_kvazaar COMMAND kvazaar_tests )
endif()

if(NOT DEFINED MSVC)
add_custom_target(
    check
    COMMAND ${PROJECT_SOURCE_DIR}/tests/test_external_symbols.sh
    COMMENT "Building and running test..."
    VERBATIM
  )
endif()