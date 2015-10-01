#ifndef STRATEGIES_COMMON_H_
#define STRATEGIES_COMMON_H_

//Use with shuffle and permutation intrinsics.
//Parameters are indices to packed elements. Each must be 0, 1, 2 or 3.
#define KVZ_PERMUTE(a, b, c, d) ( (a << 0) | (b << 2) | (c << 4) | (d << 6) )

#endif //STRATEGIES_COMMON_H_