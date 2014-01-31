#include "greatest/greatest.h"

#include "src/picture.h"


typedef struct {
  picture_list *list; // picturelist used by all tests.
} picture_list_env;

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
void picture_list_setup(picture_list_env *env)
{
  unsigned lcu_width = 2;
  unsigned lcu_height = 2;
  unsigned i;

  env->list = picture_list_init(10);
  //ASSERT(g_pl->size == 10);
  //ASSERT(g_pl->used_size == 0);

  for (i = 0; i < env->list->size; ++i) {
    picture *pic = make_pic(lcu_width, lcu_height, i);
    picture_list_add(env->list, pic);
  }
  //ASSERT(g_pl->used_size == 10);
}

/**
 * Deallocate all memory allocated by picturelist_setup.
 */
void picture_list_teardown(picture_list_env *env)
{
  picture_list_destroy(env->list);
}

//////////////////////////////////////////////////////////////////////////
// TEST FUNCTIONS

/**
 * Check that pictures have been added to the list in the correct order and
 * that they can be retrieved.
 */
TEST test_add(picture_list_env *env)
{
  unsigned i;
  for (i = 0; i < env->list->used_size; ++i) {
    picture *pic = env->list->pics[i];
    unsigned luma_size = pic->width * pic->height;
    unsigned chroma_size = luma_size >> 2;
    
    // Identify that the correct picture is in the correct place by checking
    // that the data values are the same as the index.
    ASSERT(pic->y_data[0] == i);
    ASSERT(pic->u_data[0] == i);
    ASSERT(pic->v_data[0] == i);
    ASSERT(pic->y_data[luma_size - 1] == i);
    ASSERT(pic->u_data[chroma_size - 1] == i);
    ASSERT(pic->v_data[chroma_size - 1] == i);
  }
  PASS();
}

//////////////////////////////////////////////////////////////////////////
// TEST FIXTURES
SUITE(picture_list_tests)
{
  picture_list_env env;
  SET_SETUP(picture_list_setup, &env);
  SET_TEARDOWN(picture_list_teardown, &env);

  RUN_TEST1(test_add, &env);
}
