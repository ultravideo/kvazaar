#!/usr/bin/env python3

import glob
import gzip
import os
import re
import subprocess
import sys
import tempfile
import threading
import time

# You should change these to your liking
n_threads   = 8
datadirs    = "/tmp/rdcost/data/"
resultdir   = "/tmp/rdcost/coeff_buckets"

gzargs      = ["gzip", "-d"]
filtargs    = ["./frcosts_matrix"]
octargs     = ["octave-cli", "invert_matrix.m"]
filt2args   = ["./ols_2ndpart"]

class MultiPipeManager:
    pipe_fn_template  = "%02i.txt"

    def __init__(self, odpath, dest_qps):
        self.odpath = odpath
        self.dest_qps = dest_qps

        self.pipe_fns  = []
        for qp in dest_qps:
            pipe_fn  = os.path.join(self.odpath, self.pipe_fn_template % qp)
            self.pipe_fns.append(pipe_fn)

    def __enter__(self):
        os.makedirs(self.odpath, exist_ok=True)
        for pipe_fn in self.pipe_fns:
            try:
                os.unlink(pipe_fn)
            except FileNotFoundError:
                pass
            os.mkfifo(pipe_fn)
        return self

    def __exit__(self, *_):
        for pipe_fn in self.pipe_fns:
            os.unlink(pipe_fn)

    def items(self):
        for pipe_fn in self.pipe_fns:
            yield pipe_fn

class MTSafeIterable:
    def __init__(self, iterable):
        self.lock = threading.Lock()
        self.iterable = iterable

    def __iter__(self):
        return self

    def __next__(self):
        with self.lock:
            return next(self.iterable)

def read_in_blocks(f):
    BLOCK_SZ = 65536
    while True:
        block = f.read(BLOCK_SZ)
        if (len(block) == 0):
            break
        else:
            yield block

def exhaust_gzs(sink_f, gzs):
    for gz in gzs:
        with gzip.open(gz, "rb") as f:
            if (gz == "/tmp/rdcost/data/RaceHorses_416x240_30.yuv-qp22/20.txt.gz"):
                print("kjeh")
            print("  Doing %s ..." % gz)
            for block in read_in_blocks(f):
                sink_f.write(block)
                sink_f.flush()

def run_job(jobname, input_gzs):
    resultpath = os.path.join(resultdir, "%s.result" % jobname)
    print("Running job %s" % jobname)

    with tempfile.NamedTemporaryFile() as tf:
        filt = subprocess.Popen(filtargs, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        octa = subprocess.Popen(octargs, stdin=filt.stdout, stdout=tf)

        try:
            exhaust_gzs(filt.stdin, input_gzs)
        except OSError as e:
            print("OSError %s" % e, file=sys.stderr)
            raise

        filt.stdin.close()
        filt.wait()
        octa.wait()

        if (filt.returncode != 0):
            print("First stage failed: %s" % jobname, file=sys.stderr)
            assert(0)

        with open(resultpath, "w") as rf:
            f2a = filt2args + [tf.name]
            f2 = subprocess.Popen(f2a, stdin=subprocess.PIPE, stdout=rf)
            exhaust_gzs(f2.stdin, input_gzs)
            f2.communicate()
            if (filt.returncode != 0):
                print("Second stage failed: %s" % jobname, file=sys.stderr)
                assert(0)

    print("Job %s done" % jobname)

def threadfunc(joblist):
    for jobname, job in joblist:
        run_job(jobname, job)

def scan_datadirs(path):
    seq_names = set()
    for dirent in os.scandir(path):
        if (not dirent.is_dir()):
            continue
        match = re.search("^([A-Za-z0-9_]+\.yuv)-qp[0-9]{1,2}$", dirent.name)
        if (not match is None):
            seq_name = match.groups()[0]
            seq_names.add(seq_name)

    for seq_name in seq_names:
        seq_glob = os.path.join(path, seq_name + "-qp*/")

        for qp in range(51):
            job_name = seq_name + "-qp%02i" % qp
            qp_fn = "%02i.txt.gz" % qp
            yield job_name, glob.glob(os.path.join(seq_glob, qp_fn))

def main():
    for d in (datadirs, resultdir):
        os.makedirs(d, exist_ok=True)

    jobs = scan_datadirs(datadirs)
    joblist = MTSafeIterable(iter(jobs))

    threads = [threading.Thread(target=threadfunc, args=(joblist,)) for _ in range(n_threads)]
    for thread in threads:
        thread.start()

    for thread in threads:
        thread.join()

if (__name__ == "__main__"):
    main()
