import os
import sys
import pysam
from itertools import chain, tee, izip
import gzip
import cPickle as pickle
from collections import defaultdict, Counter, namedtuple

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

def concat_files(input_list, output_path, appendnl=True):
  with open(output_path, 'w') as outf:
    for path in input_list:
      with open(path) as inf:
        nonempty = False
        for line in inf:
          nonempty = True
          outf.write(line)
        if appendnl and nonempty and not line.endswith('\n'):
          print 'WARNING last char not newline', path
          outf.write('\n')
  return

def touch(path, times=None):
  with open(path, 'a'):
    os.utime(path, times)

#--------------------------------------------------------------------------
# pickle
#--------------------------------------------------------------------------
def write_pickle(path, obj):
  f = open(path,'w')
  pickle.dump(
    obj,
    f,  
    pickle.HIGHEST_PROTOCOL
  )
  f.close()

def load_pickle(path):
  f = open(path,'r')
  obj = pickle.load(f)
  f.close()
  return obj

#--------------------------------------------------------------------------
# fa + fastq + bam
#--------------------------------------------------------------------------
def grouped(iterator, n, slop=False):
  iterator = iter(iterator)
  while True:
    vals = tuple(next(iterator, None) for _ in xrange(n))
    if None not in vals:
      yield vals
    else:
      if slop:
        vals_f = tuple(filter(lambda(x): x != None, vals))
        yield vals_f
      raise StopIteration

def fastq_iter(f):
  for lines in grouped(f, 4):
    qname, _, bcode = lines[0].strip().split('\t')[:3]
    qname = qname[1:]
    bcode = bcode.split(':')[-1]
    yield bcode, qname, lines
  raise StopIteration

def get_fasta_sizes(fa_path):
  fasta = pysam.FastaFile(fa_path)
  ctg_size_map = {}
  for ctg in fasta.references:
    size = fasta.get_reference_length(ctg)
    ctg_size_map[ctg] = size
  return ctg_size_map

def get_barcode(read):
  filt_list = filter(lambda(k, v): k == 'BC', read.tags)
  if filt_list == []: 
    return None
  else:
    k, v = filt_list[0]
    return v
  
#--------------------------------------------------------------------------
# id 
#--------------------------------------------------------------------------
class IdGenerator:
  def __init__(self, start=0):
    self.counter = start

  def get_next(self):
    i = self.counter
    self.counter += 1
    return i

  def set_counter(self, n): 
    self.counter = n 

