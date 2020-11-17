Build Kvazaar as usual with make, then edit extract_rdcosts.py so that the
parameters suit your usage (the directories, num of threads and Kvazaar
params) and then run extract_rdcosts.py. It will run a lot of Kvazaar
instances in parallel to encode a lot of videos and sift off all the coeff
groups they measure RD cost for. The coeff groups will be written into the
relevant data file in the following format (although through GZIP):

Size (B)  | Description
----------+------------
4         | size:   Coeff group size, in int16's
4         | ccc:    Coeff group's coding cost
size * 2  | coeffs: Coeff group data

You can roll your own filter_rdcosts.c program to analyze the data the way
you want, and run it like:

$ gzip -d < /path/to/compressed_datafile.gz | ./filter_rdcosts | less

Maybe one day, there'll be a multithreaded script like extract_rdcosts.py to
automate and parallelize processing of a massive heap of data files.

EDIT:
It's now possible to do OLS regression by streaming the source data twice
from source and using Octave to invert the temporary result matrix, and
that's what run_filter.py does in parallel. To do this on data you've
gathered by extract_rdcosts.py:

$ gcc filter_rdcosts.c -o frcosts_matrix
$ gcc ols_2ndpart.c -o ols_2ndpart
$ ./run_filter.py

Although you should probably adjust the run_filter.py params before actually
running it
