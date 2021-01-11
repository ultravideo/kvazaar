To extract the block costs, build Kvazaar as usual, and edit relevant
parameters in the beginning of extract_rdcosts.py and run_filter.py, most
importantly the number of cores and the set of video sequences you want to
encode to extract costs. Run extract_rdcosts.py, it will use Kvazaar to encode
each sequence and extract the costs measured there for the quantized blocks.
The costs are stored compressed and sorted by block QP, in the following
format:

Size (B)  | Description
----------+------------
4         | size:   Coeff group size, in int16's
4         | ccc:    Coeff group's coding cost
size * 2  | coeffs: Coeff group data

To analyze the costs by running a linear regression over them, build the two
tools using:

$ gcc filter_rdcosts.c -O2 -o frcosts_matrix
$ gcc ols_2ndpart.c -O2 -o ols_2ndpart

Then run the regression in parallel by running run_filter.py. The reason to do
it this way is because the data is stored compressed, so there is no way to
mmap it in Matlab/Octave/something; the data sets are absolutely huge (larger
than reasonable amounts of RAM in a decent workstation), but this way we can
store the data compressed and process it in O(1) memory complexity, so it can
be done as widely parallelized as you have CPU cores. The result files each
consist of 4 numbers, which represent an approximate linear solution to the
corresponding set of costs: the price in bits of a coefficient whose absolute
value is a) 0, b) 1, c) 2, d) 3 or higher.

After that, run rdcost_do_avg.py. It will calculate a per-QP average of the
costs over the set of the sequences having been run (ie. for each QP, take the
results for that QP for each sequence, and calculate their average). This data
is what you can use to fill in the default_fast_coeff_cost_wts table in
src/fast_coeff_cost.h.
