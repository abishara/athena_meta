import os
import sys
import cPickle as pickle
from itertools import groupby

import util

class FastqIndex(object):

  file_suffix = '.fqidx.p'

  @staticmethod
  def get_index_path(fq_path):
    return fq_path + FastqIndex.file_suffix

  @property
  def bcodes(self):
    if self._bcodes == None:
      self._bcodes = set(self._bcode_off_map.keys())
    return self._bcodes

  @property
  def num_bcodes(self):
    return len(self.bcodes)

  def __init__(
    self,
    fq_path,
    logger=None,
  ):
    
    self.logger = logger
    self.fq_path = fq_path
    self.index_path = self.get_index_path(fq_path)
    self._bcodes = None
    self._bcode_off_map = None
    self.num_pe = 0
    self.num_pe_bcoded = 0

    if not os.path.isfile(self.index_path):
      self.__build_index__()
    else:
      self.__load_index__()

    self.f_map = None
    self.open()

  def open(self):
    assert self.f_map == None, "fp map already populated"
    self.f_map = {}
    self.f_map[self.fq_path] = open(self.fq_path)
    return self

  def close(self):
    for f in self.f_map.values():
      f.close()
    return

  def __enter__(self):
    return self

  def __exit__(self, exc_type, exc_value, traceback):
    self.close()

  def __build_index__(self):  
    numbytes = 0
    self._bcode_off_map = {}
    num_pe = 0
    num_pe_bcoded = 0

    assert not self.fq_path.endswith('.gz'), \
      "gzipped fq not supported"
    with open(self.fq_path) as f:
      seen_set = set()
      for bcode, reads_iter in groupby(
        util.fastq_iter(f),
        lambda(x): x[0],
      ):
        assert bcode == None or bcode not in seen_set, \
"fastq {} NOT in barcode sorted order. Ensure reads that share barcodes \
are in a block together".format(self.fq_path)
        seen_set.add(bcode)
        if bcode != None and bcode not in self._bcode_off_map:
          self._bcode_off_map[bcode] = numbytes
        bcode_num_pe = 0
        for _, qname, lines in reads_iter:
          bcode_num_pe += 1
          txt = ''.join(lines)
          numbytes += len(txt)
        num_pe += bcode_num_pe
        if bcode != None:
          num_pe_bcoded += bcode_num_pe

    self.num_pe = num_pe
    self.num_pe_bcoded = num_pe_bcoded
    num_bcodes = len(filter(
      lambda(b): b.endswith('-1'),
      self._bcode_off_map.keys(),
    ))
    assert num_bcodes > 0, \
      "no barcodes specified in fastq {}".format(self.fq_path)
    assert self.num_pe * 0.8 < self.num_pe_bcoded, \
'''
      lower than expected ({:2.2f}%) of barcoded reads detected in fastq {}
        * use --force-input-fq to proceed
'''.format(100.*self.num_pe_bcoded / self.num_pe, self.fq_path)

    self.logger.log('fqinfo${},{},{}'.format(
      num_pe, len(self._bcode_off_map), num_bcodes,
    ))
    print 'writing index for fqs'
    for fq_path in [self.fq_path]:
      print '  -', fq_path
    util.write_pickle(self.index_path, self._bcode_off_map)

  def __load_index__(self):  
    self._bcode_off_map = util.load_pickle(self.index_path)

  def get_reads(self, bcode):
    assert self._bcode_off_map != None, 'index {} not loaded'.format(self.index_path)

    if bcode not in self.bcodes:
      raise StopIteration

    offset = self._bcode_off_map[bcode]
    f = self.f_map[self.fq_path]
    f.seek(offset, 0)

    for e in util.fastq_iter(f):
      _bcode = e[0]
      if _bcode != bcode:
        break
      yield e

    raise StopIteration

