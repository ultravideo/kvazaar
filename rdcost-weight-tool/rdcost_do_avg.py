#!/usr/bin/env python3

import glob
import sys

result_path_template = "/tmp/rdcost/coeff_buckets/*-qp%02i.result"

def main():
    results = []
    for qp in range(51):
        curr_sums = [0.0] * 4
        curr_count = 0
        result_files = glob.glob(result_path_template % qp)
        for fn in result_files:
            with open(fn) as f:
                contents = f.readlines()
                if (len(contents) != 4):
                    print("Faulty file contents at %s, skipping" % fn, file=sys.stderr)
                    continue
                nums = tuple(map(float, contents))
                if (all(n == 0.0 for n in nums)):
                    print("All-zero file %s, skipping" % fn)
                    continue

                curr_count += 1
                for i in range(len(curr_sums)):
                    curr_sums[i] += nums[i]

        if (curr_count > 0):
            curr_avgs = tuple(curr_sum / curr_count for curr_sum in curr_sums)
        else:
            curr_avgs = (0, 0, 0, 0)

        results.append(curr_avgs)
    print("\n".join(("QP %2i: " % i + ", ".join("%.6f" for _ in range(4)) % res for i, res in enumerate(results))))

if (__name__ == "__main__"):
    main()
