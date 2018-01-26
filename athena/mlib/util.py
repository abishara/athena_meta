import os
import sys
import re
import pysam
import subprocess
from itertools import chain, tee, izip
import gzip
import cPickle as pickle
from collections import defaultdict, Counter, namedtuple

#--------------------------------------------------------------------------
# prereqs
#--------------------------------------------------------------------------
def check_prereqs():
  # flye
  try:
    output = subprocess.check_output(
      'flye --version', shell=True, stderr=subprocess.STDOUT)
    vre = re.compile(r"(\d+).(\d+)")
    vr,vp = vre.search(output).groups()
    vr,vp = int(vr),int(vp)
    if None in [vr, vp]:
      print 'version cannot be parsed from --version'
      assert False
    if vr == 0 and vp == 0:
      print 'WARNING unrecognized flye version {}.{}'.format(vr,vp)
    elif vr < 2 or (vr == 2 and  vp < 3):
      print 'flye version must be 2.3+'
      assert False, "version not up to date"
  except:
    print 'flye not installed properly'
    sys.exit(1)

  # bwa
  # cannot invoke with no args without non-zero exit status

  # samtools
  try:
    output = subprocess.check_output('samtools --version', shell=True)
    vre = re.compile(r"samtools\s+(\d+).(\d+)")
    vr,vp = map(int, vre.search(output).groups())
    if None in [vr, vp]:
      print 'version cannot be parsed from --version'
      assert False
    if vr != 1 and vp < 3:
      print 'samtools version must be 1.3+'
  except:
    print 'samtools not install properly'
    sys.exit(1)

  # idba_ud
  # cannot invoke with no args without non-zero exit status

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

def _get_qname_info(line):
  def first_last(x): return (x[0], x[-1])
  elms = line.strip().split()
  qname = elms[0]
  info_map = dict(map(
    lambda(x): first_last(x.split(':')),
    elms[1:],
  ))
  return qname, info_map

def fastq_iter(f):
  for lines in grouped(f, 4):
    qname, info_map = _get_qname_info(lines[0])
    qname = qname[1:]
    assert not ('BC' in info_map and 'BX' in info_map), \
      'fastq can only either BX or BC tag!'
    bcode = None
    if 'BX' in info_map:
      bcode = info_map['BX']
    elif 'BC' in info_map:
      bcode = info_map['BC']
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
  filt_list = filter(lambda(k, v): k in ['BC', 'BX'], read.tags)
  if filt_list == []: 
    return None
  else:
    keys = set(map(lambda(k,v): k, filt_list))
    assert not ('BC' in keys and 'BX' in keys), \
      'bam can only specify either BX or BC tag!'
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

