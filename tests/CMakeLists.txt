file( GLOB TEST_SOURCES "*.c" )

# ToDo: fix the tests
list(REMOVE_ITEM TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/inter_recon_bipred_tests.c")

add_executable(kvazaar_tests ${TEST_SOURCES} )

target_include_directories(kvazaar_tests PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(kvazaar_tests PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_include_directories(kvazaar_tests PUBLIC ${PROJECT_SOURCE_DIR}/src/extras)

add_definitions(-DKVZ_DLL_EXPORTS)

if(BUILD_SHARED_LIBS)
  add_definitions(-DPIC)
endif()

list(APPEND ALLOW_AVX2 "x86_64" "AMD64")

if(MSVC)
  target_include_directories(kvazaar_tests PUBLIC ../src/threadwrapper/include)

  if(${CMAKE_SYSTEM_PROCESSOR} IN_LIST ALLOW_AVX2) 
    set_property( SOURCE ${TEST_SOURCES}  APPEND PROPERTY COMPILE_FLAGS "/arch:AVX2" )
  endif()
  add_definitions(-DWIN32_LEAN_AND_MEAN -D_WIN32 -DWIN32 -DWIN64)
else()
  if(${CMAKE_SYSTEM_PROCESSOR} IN_LIST ALLOW_AVX2) 
    set_property( SOURCE ${TEST_SOURCES}  APPEND PROPERTY COMPILE_FLAGS "-mavx2 -mbmi -mpopcnt -mlzcnt -mbmi2" )
  endif()

  if(NOT ANDROID)
    find_package(Threads REQUIRED)
    target_link_libraries(kvazaar_tests PUBLIC Threads::Threads)
  endif()

  include(CheckLibraryExists)

  CHECK_LIBRARY_EXISTS(m sin "" HAVE_LIB_M)

  if (HAVE_LIB_M)
      set(EXTRA_LIBS ${EXTRA_LIBS} m)
  endif (HAVE_LIB_M)

  target_link_libraries(kvazaar_tests PUBLIC ${EXTRA_LIBS})
endif()

if(BUILD_SHARED_LIBS)
  # The tests use internal functions, which the shared library does not
  # export, so link them to a static archive of the same objects. An archive
  # rather than the objects themselves, because some tests include library
  # sources (e.g. inter.c) and replace those objects.
  add_library(kvazaar_tests_lib STATIC)
  target_link_libraries(kvazaar_tests_lib PUBLIC kvazaar_objects)
  foreach(depth_target IN LISTS KVZ_EXTRA_BIT_DEPTH_TARGETS)
    target_sources(kvazaar_tests_lib PRIVATE $<TARGET_OBJECTS:${depth_target}>)
  endforeach()
  set_target_properties(kvazaar_tests_lib PROPERTIES LINKER_LANGUAGE C)
  target_link_libraries(kvazaar_tests PUBLIC kvazaar_tests_lib)
else()
  target_link_libraries(kvazaar_tests PUBLIC kvazaar)
endif()
