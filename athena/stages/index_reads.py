import abc
import os
import pysam
import subprocess
from collections import defaultdict
import glob

from .step import StepChunk
from ..mlib import util
from ..mlib.fq_idx import FastqIndex

class IndexReadsStep(StepChunk):

  @staticmethod
  def get_steps(options):
    # strip over fastqs to load all fq fragments
    rootfq_path = options.input_fqs
    found = False
    for fq_path in glob.glob(options.input_fqs):
      found = True
      yield IndexReadsStep(options, fq_path)
    assert found, "fastqs {} not found".format(options.input_fqs)

  def outpaths(self, final=False):
    paths = {}
    paths['pass.file'] = os.path.join(self.outdir, 'pass')
    paths['index.file'] = FastqIndex.get_index_path(self.nfq_path)
    #paths['shit.file'] = os.path.join(self.outdir, 'shit')
    return paths
 
  @property
  def outdir(self):
    return os.path.join(
      self.options.results_dir,
      self.__class__.__name__,
      str(self),
    )

  def __init__(
    self,
    options,
    fq_path,
  ):
    self.options = options
    self.fq_path = fq_path
    #self.nfq_path = fq_path[:-3]
    self.nfq_path = fq_path
    util.mkdir_p(self.outdir)

  def __fqid(self):
    return os.path.basename(os.path.dirname(os.path.dirname(self.fq_path)))

  def __str__(self):
    return '{}_{}'.format(
      self.__class__.__name__,
      self.__fqid(),
    )

  def run(self):
    #self.logger.log('uncompressing fastq')
    #cmd = 'zcat {} > {}'.format(self.fq_path, self.nfq_path)
    ##os.system(cmd)

    self.logger.log('index fastq {}'.format(self.nfq_path))

    with FastqIndex(self.nfq_path, self.logger) as idx:
      pass
    passfile_path = os.path.join(self.outdir, 'pass')
    util.touch(passfile_path)
    self.logger.log('done')

