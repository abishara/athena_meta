import os
import pysam
import random
from collections import Counter

from .step import StepChunk
from ..mlib import util
from ..mlib.fq_idx import FastqIndex

MIN_SEED_SIZE = 400
MIN_COV = 10.

class CheckReadsStep(StepChunk):

  @staticmethod
  def get_steps(options):
    assert os.path.isfile(options.input_fqs), \
      "fastqs {} not found".format(options.input_fqs)
    yield CheckReadsStep(options, options.input_fqs)

  def outpaths(self, final=False):
    paths = {}
    paths['pass.file'] = os.path.join(self.outdir, 'pass')
    paths['index.file'] = FastqIndex.get_index_path(self.nfq_path)
    paths['bins.p'] = self.options.bins_pickle_path
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
    self.logger.log('index fastq {}'.format(self.nfq_path))
    with FastqIndex(self.nfq_path, self.logger) as idx:
      fq_num_pe_bcoded = idx.num_pe_bcoded
      # check for barcodes in fastq
      assert idx.num_bcodes > 0, \
        "no barcodes specified in fastq {}".format(self.fq_path)
      assert idx.num_pe * 0.8 < idx.num_pe_bcoded, \
'''
      lower than expected ({:2.2f}%) of barcoded reads detected in fastq {}
        * use --force-reads to proceed
'''.format(100.*idx.num_pe_bcoded / idx.num_pe, self.fq_path)

    # use cheat seeds if specified (for debugging)
    if self.options.cheat_seeds:
      self.logger.log('using cheat seeds file: {}'.format(self.options.cheat_seeds))
      seeds = set()
      with open(self.options.cheat_seeds) as fin:
        for line in fin:
          seed = line.strip()
          seeds.add(seed)
      self.logger.log('  - loaded {}'.format(len(seeds)))
      seeds = list(seeds)
    # use read mappings from *bam to select seeds without high enough input
    # coverage
    else:
      self.logger.log('get seed contigs from input assembly')
      ctg_covs, bam_num_se_bcoded = self.get_bam_stats()
      bam_num_pe_bcoded = bam_num_se_bcoded / 2
      assert bam_num_pe_bcoded < 0.8 * fq_num_pe_bcoded, \
'''
      only ~{:2.2f}% of barcoded reads detected in *bam {} as compared to
      fastq {}
        * use --force-reads to proceed
'''.format(100.*bam_num_pe_bcoded / fq_num_pe_bcoded, self.fq_path)

      seeds = self.get_seeds(ctg_covs)
      random.shuffle(seeds)

    # strip seed contigs into bins such that no more than 4000 bins
    bins = []
    group_size = max(1, len(seeds) / 4000)
    for i, seed_group in \
      enumerate(util.grouped(seeds, group_size, slop=True)):
      binid = 'bin.{}'.format(i)
      bins.append((binid, seed_group))

    self.logger.log('created {} bins from seeds'.format(len(bins)))
    util.write_pickle(self.options.bins_pickle_path, bins)

    passfile_path = os.path.join(self.outdir, 'pass')
    util.touch(passfile_path)
    self.logger.log('done')

  def get_seeds(self, ctg_covs):
    ctg_size_map = util.get_fasta_sizes(self.options.ctgfasta_path)
    seeds = ctg_size_map.keys()
    seeds = filter(
      lambda(c): (
        ctg_size_map[c] >= MIN_SEED_SIZE and
        ctg_covs[c] >= MIN_COV
      ),
      seeds,
    )
    self.logger.log('  {} total inputs seeds covering {} bases'.format(
      len(ctg_size_map), sum(ctg_size_map.values())
    ))
    self.logger.log('  {} input seed contigs >= {}bp and >= {}x coverage covering {} bases'.format(
      len(seeds),
      MIN_SEED_SIZE,
      MIN_COV,
      sum(map(lambda(c): ctg_size_map[c], seeds)),
    ))
    return seeds

  def get_bam_stats(self):
    ctg_counts_path = os.path.join(self.options.working_dir, 'ctg_counts.p')
    if os.path.isfile(ctg_counts_path):
      return Counter(util.load_pickle(ctg_counts_path))
    self.logger.log('computing seed coverages (required pass thru *bam)')
    bam_fin = pysam.Samfile(self.options.reads_ctg_bam_path, 'rb')
    ctg_bases = Counter()
    
    num_se_bcoded = 0
    for i, read in enumerate(bam_fin):
      if not read.is_secondary and util.get_barcode(read) != None:
        num_se_bcoded += 1
      if read.is_unmapped:
        continue
      seed_ctg = bam_fin.getrname(read.tid) 
      ctg_bases[seed_ctg] += read.query_alignment_length
    ctg_size_map = util.get_fasta_sizes(self.options.ctgfasta_path)
    ctg_covs = Counter(dict(map(
      lambda(c, b) : (c, 1. * b / ctg_size_map[c]),
      ctg_bases.iteritems()
    )))
    util.write_pickle(ctg_counts_path, ctg_covs)
    return ctg_covs, num_se_bcoded

