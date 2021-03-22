#include "fast_coeff_cost.h"
#include "kvazaar.h"

// Note: Assumes that costs are non-negative, for pretty obvious reasons
static uint16_t to_q88(float f)
{
  return (uint16_t)(f * 256.0f + 0.5f);
}

static uint64_t to_4xq88(const float f[4])
{
  int i;
  uint64_t result = 0;

  for (i = 3; i >= 0; i--) {
    result <<= 16;
    result |= to_q88(f[i]);
  }
  return result;
}

int kvz_fast_coeff_table_parse(fast_coeff_table_t *fast_coeff_table, FILE *fast_coeff_table_f)
{
  int i;
  uint64_t *wts_by_qp = fast_coeff_table->wts_by_qp;

  for (i = 0; i < MAX_FAST_COEFF_COST_QP; i++) {
    float curr_wts[4];

    if (fscanf(fast_coeff_table_f, "%f %f %f %f\n", curr_wts + 0,
                                                    curr_wts + 1,
                                                    curr_wts + 2,
                                                    curr_wts + 3) != 4) {
      return 1;
    }
    wts_by_qp[i] = to_4xq88(curr_wts);
  }
  return 0;
}

void kvz_fast_coeff_use_default_table(fast_coeff_table_t *fast_coeff_table)
{
  int i;
  uint64_t *wts_by_qp = fast_coeff_table->wts_by_qp;

  for (i = 0; i < MAX_FAST_COEFF_COST_QP; i++) {
    wts_by_qp[i] = to_4xq88(default_fast_coeff_cost_wts[i]);
  }
}

uint64_t kvz_fast_coeff_get_weights(const encoder_state_t *state)
{
  const fast_coeff_table_t *table = &(state->encoder_control->fast_coeff_table);
  return table->wts_by_qp[state->qp];
}
