#!/usr/bin/env python3

import glob
import gzip
import os
import subprocess
import threading
import time

logdir = os.path.join("/tmp", "rdcost", "logs")
ofdir  = os.path.join("/tmp", "rdcost", "data")

# Note that n_kvazaars * len(dest_qps) has to be less than the max number of
# fd's that a process can have (check it out: ulimit -a, likely 1024)
smt_threads   = 8 # Kinda lazy, but just match this to your cpu
n_kvz_threads = 1 # How many threads each kvz instance is running?
n_kvazaars    = smt_threads // n_kvz_threads
kvz_srcdir    = lambda path: os.path.join(os.path.dirname(os.path.realpath(__file__)), path)

dest_qps      = tuple(range(51))
base_qps      = tuple(range(22, 28))
sequences     = ("/opt/test_seqs/hevc-D/*.yuv",)# "/opt/test_seqs/custom_seqs/*/*.yuv")

kvzargs       = [kvz_srcdir("kvazaar"), "--threads", str(n_kvz_threads), "--preset=ultrafast", "--fastrd-sampling", "--fast-residual-cost=0"]
kvzenv        = {"LD_LIBRARY_PATH": kvz_srcdir(".libs/")}

class MultiPipeGZOutManager:
    pipe_fn_template  = "%02i.txt"
    gzout_fn_template = "%02i.txt.gz"

    def __init__(self, odpath, dest_qps):
        self.odpath = odpath
        self.dest_qps = dest_qps

        self.pipe_fns  = []
        self.gzout_fns = []
        for qp in dest_qps:
            pipe_fn  = os.path.join(self.odpath, self.pipe_fn_template % qp)
            gzout_fn = os.path.join(self.odpath, self.gzout_fn_template % qp)

            self.pipe_fns.append(pipe_fn)
            self.gzout_fns.append(gzout_fn)

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
        for pipe_fn, gzout_fn in zip(self.pipe_fns, self.gzout_fns):
            yield (pipe_fn, gzout_fn)

class MTSafeIterable:
    def __init__(self, iterable):
        self.lock = threading.Lock()
        self.iterable = iterable

    def __iter__(self):
        return self

    def __next__(self):
        with self.lock:
            return next(self.iterable)

def combinations(xi, yi):
    for x in xi:
        for y in yi:
            yield (x, y)

def chain(lol):
    for l in lol:
        for i in l:
            yield i

# Would've used Popen with gzip, but "gzip [fifo]" with an unconnected fifo
# will detect the situation and not block, but just consider it an empty
# file. Don't like it when tools outsmart their user..
def do_gzip(in_fn, out_fn):
    BLOCK_SZ = 65536
    PRINT_MULT = 1024
    with open(in_fn, "rb") as inf, gzip.open(out_fn, "wb") as outf:
        num_read = 0
        print_next_thres = BLOCK_SZ * PRINT_MULT
        while True:
            block = inf.read(BLOCK_SZ)
            num_read += len(block)
            if (num_read >= print_next_thres):
                print("    read     %8i MB from %s" % (num_read / (1024 * 1024), in_fn))
                print_next_thres += BLOCK_SZ * PRINT_MULT

            if (len(block) == 0):
                break
            outf.write(block)

        print("    finished %8i MB from %s" % (num_read / (1024 * 1024), in_fn))

def run_job(job):
    ifpath, qp = job
    ifname = os.path.basename(ifpath)

    jobname  = "%s-qp%i" % (ifname, qp)
    hevcname = "%s.hevc" % jobname
    logname  = "%s.log"  % jobname
    odname   = jobname

    hevcpath = os.path.join("/tmp", hevcname)
    logpath  = os.path.join(logdir, logname)
    odpath   = os.path.join(ofdir,  odname)

    my_kvzargs = kvzargs + ["-i",              ifpath,
                            "--qp",            str(qp),
                            "-o",              hevcpath,
                            "--fastrd-outdir", odpath]

    with open(logpath, "w") as lf:
        with MultiPipeGZOutManager(odpath, dest_qps) as pipes_and_outputs:
            gzips = []
            gzip_threads = []
            for pipe_fn, out_fn in pipes_and_outputs.items():
                gzip_thread = threading.Thread(target=do_gzip, args=(pipe_fn, out_fn))
                gzip_thread.start()
                gzip_threads.append(gzip_thread)

            kvz = subprocess.Popen(my_kvzargs, env=kvzenv, stderr=lf)

            kvz.communicate()
            for gzip in gzips:
                gzip.communicate()

def threadfunc(joblist):
    for job in joblist:
        run_job(job)

def main():
    jobs = combinations(chain(map(glob.glob, sequences)), base_qps)
    joblist = MTSafeIterable(jobs)

    threads = [threading.Thread(target=threadfunc, args=(joblist,)) for _ in range(n_kvazaars)]
    for thread in threads:
        thread.start()

    for thread in threads:
        thread.join()

if (__name__ == "__main__"):
    main()
