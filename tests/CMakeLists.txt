file( GLOB TEST_SOURCES "*.c" )

# ToDo: fix the tests
list(REMOVE_ITEM TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/inter_recon_bipred_tests.c")

add_executable(kvazaar_tests ${TEST_SOURCES} )

target_include_directories(kvazaar_tests PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(kvazaar_tests PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_include_directories(kvazaar_tests PUBLIC ${PROJECT_SOURCE_DIR}/src/extras)

if(USE_SHARED_LIB)
  message(INFO " tests do not work with shared lib at the moment")
  add_definitions(-DPIC -DKVZ_DLL_EXPORTS)
endif()

if(MSVC)
  target_include_directories(kvazaar_tests PUBLIC ../src/threadwrapper/include)

  set_property( SOURCE ${TEST_SOURCES}  APPEND PROPERTY COMPILE_FLAGS "/arch:AVX2" )
  add_definitions(-DWIN32_LEAN_AND_MEAN -D_WIN32 -DWIN32 -DWIN64)
else()
  set_property( SOURCE ${TEST_SOURCES}  APPEND PROPERTY COMPILE_FLAGS "-mavx2 -mbmi -mpopcnt -mlzcnt -mbmi2" )
  find_package(Threads REQUIRED)
  target_link_libraries(kvazaar_tests PUBLIC Threads::Threads)

  include(CheckLibraryExists)

  CHECK_LIBRARY_EXISTS(m sin "" HAVE_LIB_M)

  if (HAVE_LIB_M)
      set(EXTRA_LIBS ${EXTRA_LIBS} m)
  endif (HAVE_LIB_M)

  target_link_libraries(kvazaar_tests PUBLIC ${EXTRA_LIBS})
endif()

add_definitions(-DCOMPILE_INTEL)

target_link_libraries(kvazaar_tests PUBLIC libuvg266)

