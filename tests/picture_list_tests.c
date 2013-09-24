#include "seatest.h"

#include "picture.h"

//////////////////////////////////////////////////////////////////////////
// GLOBALS
picture_list *g_pl; // picturelist used by all tests.


//////////////////////////////////////////////////////////////////////////
// SETUP, TEARDOWN AND HELPER FUNCTIONS

/**
 * Alloc and initialize picture. Initialize yuv-data with parameter init.
 */
picture * make_pic(unsigned lcu_width, unsigned lcu_height, char init) {
  picture *pic = picture_init(LCU_WIDTH * lcu_width, LCU_WIDTH * lcu_height, 
    lcu_width, lcu_height);
  unsigned luma_size = pic->width * pic->height;
  unsigned chroma_size = luma_size >> 2;

  memset(pic->y_data, init, luma_size);
  memset(pic->u_data, init, chroma_size);
  memset(pic->v_data, init, chroma_size);

  return pic;
}

/**
 * Setup initial conditions for all tests.
 * 
 * The conditions are that g_pl is a picture list of size 10 with 10 pictures.
 * The pictures have their yuv-data initialized with their index in the list,
 * meaning that the first picture is filled with 0 values and the second with
 * 1 values.
 */
void picturelist_setup(void)
{
  unsigned lcu_width = 2;
  unsigned lcu_height = 2;
  unsigned i;

  g_pl = picture_list_init(10);
  assert_true(g_pl->size == 10);
  assert_true(g_pl->used_size == 0);

  for (i = 0; i < g_pl->size; ++i) {
    picture *pic = make_pic(lcu_width, lcu_height, i);
    picture_list_add(g_pl, pic);
  }
  assert_true(g_pl->used_size == 10);
}

/**
 * Deallocate all memory allocated by picturelist_setup.
 */
void picturelist_teardown(void)
{
  picture_list_destroy(g_pl);
}

//////////////////////////////////////////////////////////////////////////
// TEST FUNCTIONS

/**
 * Check that pictures have been added to the list in the correct order and
 * that they can be retrieved.
 */
void test_add(void)
{
  unsigned i;
  for (i = 0; i < g_pl->used_size; ++i) {
    picture *pic = g_pl->pics[i];
    unsigned luma_size = pic->width * pic->height;
    unsigned chroma_size = luma_size >> 2;
    
    // Identify that the correct picture is in the correct place by checking
    // that the data values are the same as the index.
    assert_true(pic->y_data[0] == i);
    assert_true(pic->u_data[0] == i);
    assert_true(pic->v_data[0] == i);
    assert_true(pic->y_data[luma_size - 1] == i);
    assert_true(pic->u_data[chroma_size - 1] == i);
    assert_true(pic->v_data[chroma_size - 1] == i);
  }
}

//////////////////////////////////////////////////////////////////////////
// TEST FIXTURES
void picture_list_tests(void)
{
  test_fixture_start();
  fixture_setup(picturelist_setup);
  run_test(test_add);
  fixture_teardown(picturelist_teardown);
  test_fixture_end();
}
