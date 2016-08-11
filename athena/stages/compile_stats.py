import os
import pysam
import subprocess
from collections import Counter

from .step import StepChunk
from ..mlib import util

class CompileStatsStep(StepChunk):

  @classmethod
  def get_steps(cls, options):

    ctg_size_map = util.get_fasta_sizes(options.ctgfasta_path)
    if os.path.isfile(options.bins_pickle_path):
      bins = util.load_pickle(options.bins_pickle_path)
    else:
      bins = map(
        lambda(b, v): ('bin.{}'.format(b), v),
        enumerate(util.get_fasta_partitions(options.ctgfasta_path, 128)),
      )
      print '  - saving {}'.format(options.bins_pickle_path)
      util.write_pickle(options.bins_pickle_path, bins)

    total_ctgs = 0
    total_bases = 0
    for binid, group in bins:
      total_bases += sum(map(lambda(c): ctg_size_map[c], group))
      total_ctgs += len(group)
      yield cls(options, binid, group)

    print 'grouped {} total_ctgs covering {} bases'.format(
      total_ctgs,
      total_bases)

  @property
  def outdir(self):
    return os.path.join(
      self.options.results_dir,
      self.__class__.__name__,
      str(self),
    )

  def outpaths(self, final=False):
    paths = {}
    paths['stats.p'] = os.path.join(self.outdir, 'stats.p')
    return paths
 
  def __str__(self):
    return '{}_{}'.format(
      self.__class__.__name__,
      self.binid,
    )

  def __init__(
    self,
    options,
    binid,
    ctgs,
  ):
    self.options = options
    self.ctgs = ctgs
    self.binid = binid
    util.mkdir_p(self.outdir)
    util.mkdir_p(self.options.get_bin_dir(self.binid))

  def run(self):
    # count reads for each contig in the group
    fhandle = pysam.Samfile(self.options.reads_ctg_bam_path, 'rb')
    ctg_counts = Counter()
    for ctg in self.ctgs:
      self.logger.log('counting reads on ctg {}'.format(ctg))
      for read in fhandle.fetch(ctg):
        if (
          not read.is_unmapped and 
          not read.is_secondary and 
          not read.is_supplementary
        ):
          ctg_counts[ctg] += 1
    fhandle.close()

    self.logger.log('done')
    out_path = os.path.join(self.outdir, 'stats.p')
    util.write_pickle(out_path, ctg_counts)

