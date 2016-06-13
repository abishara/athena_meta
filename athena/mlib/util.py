import os
import sys
import pysam
from itertools import chain, tee, izip
import gzip
import cPickle as pickle
from bx.intervals.intersection import IntervalTree
from collections import defaultdict, Counter

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
# intervals
#--------------------------------------------------------------------------
def load_bed(bed_path):
  intervals_map = defaultdict(lambda: IntervalTree(1,1))
  with open(bed_path) as f:
    for i, line in enumerate(f):
      if line.startswith('#'):
        continue
      words = line.split('\t')
      ctg = words[0]
      b = int(words[1])
      e = int(words[2])
      intervals_map[ctg].insert(b, e, i)
  return intervals_map

#--------------------------------------------------------------------------
# partitioning
#--------------------------------------------------------------------------
def pairwise(iterable):
  a, b = tee(iterable)
  next(b, None)
  return izip(a, b)

def get_partitions(begin, end, step):
  return pairwise(chain(xrange(begin,end,step), [end]))

def get_genome_partitions(
  fasta_path,
  window_size,
  step_size,
  sel_intervals=None,
  sel_ctgs=None,
):
  f = pysam.FastaFile(fasta_path)
  
  for ctg in f.references:
    if sel_ctgs and ctg not in sel_ctgs:
      continue
    l = f.get_reference_length(ctg)
    for (b, _) in get_partitions(0, l, step_size):
      e = b + window_size
      yield (ctg, b, e)
  f.close()
  
#--------------------------------------------------------------------------
# fastq
#--------------------------------------------------------------------------
def grouped(iterator, lines_per_read):
  while True:
    vals = tuple(next(iterator, None) for _ in xrange(lines_per_read))
    if None not in vals:
      yield vals
    else:
      raise StopIteration

def fa_iter(fa):
  with open(fa) as f:
    for fields in grouped(f, 2):
      qname_ln = fields[0]
      qname = qname_ln[1:].strip()
      txt = ''.join(fields)
      yield (qname, txt)
  raise StopIteration

def tenx_fastq_iter(fq, fmt='raw'):
  assert fmt in [
    'raw',
    'fa',
    'fa-bcode',
  ]
  assert os.path.isfile(fq)
  open_f = open
  if fq.endswith('.gz'):
    open_f = gzip.open
  with open_f(fq) as f:

    for fields in grouped(f, 9):
      bcode = fields[5].strip()
      if fmt == 'raw':
        rtxt = ''.join(fields)
        yield (bcode, rtxt)
      elif fmt in ['fa', 'fa-bcode']:
        qname_ln = fields[0]
        # strip aux fields
        qname_w = qname_ln.split()[0]
        # strip /1,/2
        qname = qname_w.split('/')[0]
        # strip @
        assert qname.startswith('@')
        qname = qname[1:]

        if fmt == 'fa-bcode':
          header = '>'+bcode+'$'+qname+'\n'
        else:
          header = '>'+qname+'\n'
        rtxt = ''.join([
          header,
          fields[1],
          header,
          fields[3],
        ])
        yield (bcode, rtxt)

  raise StopIteration

def convert_tenx_fq2fa(tenxfq_path, fa_path):
  with open(fa_path, 'w') as f:
    for _, txt in tenx_fastq_iter(tenxfq_path, fmt='fa'):
      f.write(txt)
def convert_tenx_fq2bcodefa(tenxfq_path, fa_path):
  with open(fa_path, 'w') as f:
    for _, txt in tenx_fastq_iter(tenxfq_path, fmt='fa-bcode'):
      f.write(txt)

