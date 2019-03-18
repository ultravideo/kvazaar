#!/usr/bin/env python3

import glob
import os
import subprocess
import tempfile
import threading
import time

n_threads   = 8
data        = "/home/moptim/rdcost/data/*.gz"
gzargs      = ["gzip", "-d"]
filtargs    = ["./frcosts_matrix"]
octargs     = ["octave-cli", "invert_matrix.m"]
filt2args   = ["./ols_2ndpart"]
resultdir   = os.path.join("/tmp", "rdcost", "coeff_buckets")

class MTSafeIterable:
    def __init__(self, iterable):
        self.lock = threading.Lock()
        self.iterable = iterable

    def __iter__(self):
        return self

    def __next__(self):
        with self.lock:
            return next(self.iterable)

def run_job(job):
    datapath   = job
    resultpath = os.path.join(resultdir, os.path.basename(job) + ".result")

    print("Running job %s" % datapath)

    with open(resultpath, "w") as rf:
        with tempfile.NamedTemporaryFile() as tf:
            with open(datapath, "rb") as df:
                f2a = list(filt2args)
                f2a.append(tf.name)
                gzip = subprocess.Popen(gzargs, stdin=df, stdout=subprocess.PIPE)
                filt = subprocess.Popen(filtargs, stdin=gzip.stdout, stdout=subprocess.PIPE)
                octa = subprocess.Popen(octargs, stdin=filt.stdout, stdout=tf)

                octa.communicate()
                filt.communicate()
                gzip.communicate()

            with open(datapath, "rb") as df:
                gz2  = subprocess.Popen(gzargs, stdin=df, stdout=subprocess.PIPE)
                f2   = subprocess.Popen(f2a, stdin=gz2.stdout, stdout=rf)

                f2.communicate()
                gz2.communicate()

    print("Job %s done" % datapath)

def threadfunc(joblist):
    for job in joblist:
        run_job(job)

def main():
    jobs = glob.glob(data)
    joblist = MTSafeIterable(iter(jobs))

    threads = [threading.Thread(target=threadfunc, args=(joblist,)) for _ in range(n_threads)]
    for thread in threads:
        thread.start()

    for thread in threads:
        thread.join()

if (__name__ == "__main__"):
    main()
