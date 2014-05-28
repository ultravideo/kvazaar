import numpy as np
import matplotlib.pyplot as plt

import re, weakref

class LogJob:
  def __init__(self, worker_id, enqueue, start, stop, dequeue, description, is_thread_job=True):
    self._worker_id = worker_id
    self._enqueue = enqueue
    self._start = start
    self._stop = stop
    self._dequeue = dequeue
    self._description = description
    self._is_thread_job = is_thread_job

    self._depends = []
    self._rdepends = []

  def add_dependency_on(self, o):
    self._depends.append(weakref.proxy(o))
    o._rdepends.append(weakref.proxy(self))

  def remove_offset(self, offset):
    if self._enqueue is not None:
      self._enqueue -= offset
    self._start -= offset
    self._stop -= offset
    if self._dequeue is not None:
      self._dequeue -= offset

  def get_min_offset(self):
    return min([x for x in [self._enqueue, self._start, self._stop, self._dequeue] if x is not None])

  def _get_properties(self):
    return dict([x.split('=',1) for x in self._description.split(',')])

  def _position_from_str(self, s):
    if re.match('^[0-9]+$', s):
      return int(s)
    else:
      v = [float(x) for x in s.split('-', 1)]
      return (v[0] + v[1]) / 2

  def _height_from_str(self, s):
    if self._is_thread_job:
      diff = 0.2
    else:
      diff = 0.4
    if re.match('^[0-9]+$', s):
      return 1.0 - diff
    else:
      v = [float(x) for x in s.split('-', 1)]
      return (max(v) - min(v) + 1) - diff

  def height(self):
    desc = self._get_properties()
    if 'row' in desc:
      return self._height_from_str(desc['row'])
    elif 'position_y' in desc:
      return self._height_from_str(desc['position_y'])
    else:
      if self._is_thread_job:
        return 0.8
      else:
        return 0.6

  def position_y(self):
    desc = self._get_properties()
    if 'row' in desc:
      return self._position_from_str(desc['row'])
    elif 'position_y' in desc:
      return self._position_from_str(desc['position_y'])
    else:
      return -1

class LogFlush:
  def __init__(self, when):
    self._when = when

  def remove_offset(self, offset):
    self._when -= offset

  def get_min_offset(self):
    return self._when

class LogThread:
  def __init__(self, worker_id, start, stop):
    self._worker_id = worker_id
    self._start = start
    self._stop = stop

  def remove_offset(self, offset):
    self._start -= offset
    self._stop -= offset

  def get_min_offset(self):
    return min([self._start, self._stop])

  def plot(self, ax, i):
    ax.barh(i, self._stop - self._start, left=self._start, height=0.9, align='center',label="test", color='yellow')

class LogParser:
  def _parse_time(self, base, sign, value):
    if sign == '+':
      return base + float(value)
    else:
      return float(value)

  def __init__(self, filename):
    re_thread = re.compile(r'^\t([0-9]+)\t-\t([0-9\.]+)\t(\+?)([0-9\.]+)\t-\tthread$')
    re_job = re.compile(r'^([^\t]+)\t([0-9]+)\t([0-9\.]+)\t(\+?)([0-9\.]+)\t(\+?)([0-9\.]+)\t(\+?)([0-9\.]+)\t(.*)$')
    re_dep = re.compile(r'^(.+)->(.+)$')
    re_flush = re.compile(r'^\t\t-\t-\t([0-9\.]+)\t-\tFLUSH$')
    re_other_perf = re.compile(r'^\t([0-9]*)\t-\t([0-9\.]+)\t(\+?)([0-9\.]*)\t-\t(.*)$')

    objects_cache = {}
    deps_cache = []
    objects = []
    threads = {}

    for line in open(filename,'r').readlines():
      m = re_thread.match(line)
      if m:
        g = m.groups()

        thread_id = int(g[0])
        start = self._parse_time(0, '', g[1])
        stop = self._parse_time(start, g[2], g[3])

        threads[thread_id] = LogThread(thread_id, start, stop)
        continue

      m = re_flush.match(line)
      if m:
        g = m.groups()
        when = self._parse_time(0, '', g[0])
        objects.append(LogFlush(when))

        for g in deps_cache:
          objects_cache[g[1]].add_dependency_on(objects_cache[g[0]])

        #clear object cache
        objects_cache = {}
        deps_cache = []
        continue

      m = re_job.match(line)
      if m:
        g = m.groups()
        worker_id = int(g[1])
        enqueue = self._parse_time(0, '', g[2])
        start = self._parse_time(enqueue, g[3], g[4])
        stop = self._parse_time(start, g[5], g[6])
        dequeue = self._parse_time(stop, g[7], g[8])
        description = g[9]

        value = LogJob(worker_id, enqueue, start, stop, dequeue, description)
        objects.append(value)
        objects_cache[g[0]] = value
        continue

      m = re_dep.match(line)
      if m:
        g = m.groups()
        deps_cache.append((g[0],g[1]))
        continue

      m = re_other_perf.match(line)
      if m:
        g = m.groups()
        if g[0] == '':
          worker_id = None
        else:
          worker_id = int(g[0])

        start = self._parse_time(0, '', g[1])
        stop = self._parse_time(start, g[2], g[3])

        objects.append(LogJob(worker_id, None, start, stop, None, g[4], False))
        continue

      raise ValueError("Unknown line:", line)

    assert len(threads) + len(objects) > 0
    self._threads = threads
    self._objects = objects

    #Remove offsets
    offset = min([x.get_min_offset() for x in self._threads.values()] + [x.get_min_offset() for x in self._objects])

    for x in self._threads.values() + self._objects:
      x.remove_offset(offset)

  def plot_threads(self):
    fig = plt.figure()
    ax=fig.gca()
    yticks = {}
    for k in sorted(self._threads.keys()):
      v = self._threads[k]
      v.plot(ax, -k)
      yticks[-k] = 'Thread {0}'.format(k)

    for o in self._objects:
      if isinstance(o, LogJob):
        ax.barh(-o._worker_id, o._stop - o._start, left=o._start, height=0.8, align='center',label="test", color='green')

      if isinstance(o, LogFlush):
        ax.axvline(o._when)

    for o in self._objects:
      if isinstance(o, LogJob):
        for o2 in o._depends:
          ax.plot([o2._stop, o._start], [-o2._worker_id, -o._worker_id], linewidth=2, color='r')

    plt.yticks( yticks.keys(), yticks.values() )
    fig.show()
    plt.show()

  def get_color(self, i, is_thread_job=True):
    if i is None:
      return 'w'
    if is_thread_job:
      color_keys = ['#ff0000', '#00ff00', '#0000ff', '#ffff00', '#ff00ff', '#00ffff']
    else:
      color_keys = ['#ffaaaa', '#aaffaa', '#aaaaff', '#ffffaa', '#ffaaff', '#aaffff']

    return color_keys[i%len(color_keys)]

  def plot_picture_wise_wpp(self):
    fig = plt.figure()
    ax=fig.gca()

    yticks = {}

    #first draw threads
    for o in self._objects:
      if isinstance(o, LogJob) and o._is_thread_job:
        y = o.position_y()
        height = o.height()
        ax.barh(-y, o._stop - o._start, left=o._start, height=height, align='center', color=self.get_color(o._worker_id))
        if y % 1 <= 0.0001:
          yticks[int(-y)]= int(y)

    #then jobs
    for o in self._objects:
      if isinstance(o, LogJob) and not o._is_thread_job:
        y = o.position_y()
        height = o.height()
        ax.barh(-y, o._stop - o._start, left=o._start, height=height, align='center', color=self.get_color(o._worker_id, False))
        if y % 1 <= 0.0001:
          yticks[int(-y)]= int(y)

    for o in self._objects:
      if isinstance(o, LogJob):
        for o2 in o._depends:
          ax.plot([o2._stop, o._start], [-int(o2.position_y()), -int(o.position_y())] , linewidth=1, color='k')

    for y in yticks.keys():
      ax.axhline(y+0.5, color='k')
      if y - 1 not in yticks.keys():
        ax.axhline(y-0.5, color='k')
      if y == 1:
        yticks[y] = "None"

    plt.yticks( yticks.keys(), yticks.values() )

    for o in self._objects:
      if isinstance(o, LogFlush):
        ax.axvline(o._when)


    ax.set_xlabel("Time [s]")
    ax.set_ylabel("LCU y coordinate")

    fig.show()
    plt.show()


  def plot_animation(self):
    pass


if __name__ == '__main__':
  import sys
  if len(sys.argv) > 1:
    l = LogParser(sys.argv[1])
  else:
    l = LogParser('threadqueue.log')
  l.plot_picture_wise_wpp()
  #l.plot_threads()
