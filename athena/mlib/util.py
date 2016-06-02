import os
import sys
from itertools import chain, tee, izip

#--------------------------------------------------------------------------
# os
#--------------------------------------------------------------------------
def mkdir_p(dir):
  if not os.path.isdir(dir):
    os.makedirs(dir)

def mktmpdir(prefix='tmp'):
  path = tempfile.mkdtemp(dir='.', prefix=prefix)
  return path

class cd: 
  def __init__(self, newPath):
    self.newPath = newPath

  def __enter__(self):
    self.savedPath = os.getcwd()
    os.chdir(self.newPath)

  def __exit__(self, etype, value, traceback):
    os.chdir(self.savedPath)

def concatFiles(input_list, output_path):
  with open(output_path, 'w') as outf:
    for path in input_list:
      with open(path) as inf:
        for line in inf:
          outf.write(line)

  return

def touch(fname, times=None):
  with open(fname, 'a'):
    os.utime(fname, times)

#--------------------------------------------------------------------------
# partitioning
#--------------------------------------------------------------------------
def pairwise(iterable):
  a, b = tee(iterable)
  next(b, None)
  return izip(a, b)

def get_partitions(begin, end, step):
  return pairwise(chain(xrange(begin,end,step), [end]))

