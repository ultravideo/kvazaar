#!/usr/bin/env python3

import glob
import os
import subprocess
import threading
import time

logdir = os.path.join("/tmp", "rdcost", "logs")
ofdir  = os.path.join("/tmp", "rdcost", "data")

n_threads   = 8
home_rel    = lambda path: os.path.join(os.environ["HOME"], path)

qps         = range(12, 14)
sequences   = ("/opt/test_seqs/hevc-B/*.yuv",)

kvzargs     = [home_rel("kvazaar/src/kvazaar"), "--threads", "1", "--preset=ultrafast"]
kvzenv      = {"LD_LIBRARY_PATH": home_rel("kvazaar/src/.libs/")}

gzargs      = ["gzip"]

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

def run_job(job):
    ifpath, qp = job
    ifname = os.path.basename(ifpath)

    jobname  = "%s-qp%i" % (ifname, qp)
    hevcname = "%s.hevc" % jobname
    logname  = "%s.log"  % jobname
    ofname   = "%s.gz"   % jobname

    hevcpath = os.path.join("/tmp", hevcname)
    logpath  = os.path.join(logdir, logname)
    ofpath   = os.path.join(ofdir,  ofname)

    my_kvzargs = kvzargs + ["-i",   ifpath,
                            "--qp", str(qp),
                            "-o",   hevcpath]

    with open(logpath, "w") as lf:
        with open(ofpath, "wb") as of:
            kvz = subprocess.Popen(my_kvzargs, env=kvzenv, stdout=subprocess.PIPE, stderr=lf)
            gzip = subprocess.Popen(gzargs, stdin=kvz.stdout, stdout=of)

            gzip.communicate()
            kvz.communicate()

def threadfunc(joblist):
    for job in joblist:
        run_job(job)

def main():
    jobs = combinations(chain(map(glob.glob, sequences)), qps)
    joblist = MTSafeIterable(jobs)

    threads = [threading.Thread(target=threadfunc, args=(joblist,)) for _ in range(n_threads)]
    for thread in threads:
        thread.start()

    for thread in threads:
        thread.join()

if (__name__ == "__main__"):
    main()
