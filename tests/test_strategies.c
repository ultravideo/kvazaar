#include "test_strategies.h"

#include "src/strategyselector.h"


strategy_list strategies;


void init_test_strategies()
{
  strategies.allocated = 0;
  strategies.count = 0;
  strategies.strategies = NULL;

  // Init strategyselector because it sets hardware flags.
  strategyselector_init(1);

  // Collect all strategies to be tested.
  if (!strategy_register_picture(&strategies)) {
    fprintf(stderr, "strategy_register_picture failed!\n");
    return;
  }

  if (!strategy_register_dct(&strategies)) {
    fprintf(stderr, "strategy_register_partial_butterfly failed!\n");
    return;
  }
}
