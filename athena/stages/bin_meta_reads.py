import os
import pysam
import shutil
from bx.intervals.cluster import ClusterTree
from collections import defaultdict, Counter
import numpy as np
from itertools import groupby, combinations
import random

from .step import StepChunk
from ..mlib import util

MIN_SEED_SIZE = 500
MIN_SEED_SIZE = 400
MIN_COV = 10.

class BinMetaReadsStep(StepChunk):

  @staticmethod
  def get_steps(options):
    yield BinMetaReadsStep(options)

  def __init__(
    self,
    options,
  ):
    self.options = options
    util.mkdir_p(self.outdir)

  def __str__(self):
    return self.__class__.__name__

  @property
  def outdir(self):
    return os.path.dirname(self.options.groups_pickle_path)

  def outpaths(self, final=False):
    paths = {}
    paths['bins.p'] = self.options.bins_pickle_path
    #paths['shit.p'] = 'shit'
    return paths

  def get_ctg_read_counts(self):
    ctg_counts_path = os.path.join(self.options.working_dir, 'ctg_counts.p')
    if os.path.isfile(ctg_counts_path):
      return Counter(util.load_pickle(ctg_counts_path))
    self.logger.log('computing seed coverages (required pass thru *bam)')
    bam_fin = pysam.Samfile(self.options.reads_ctg_bam_path, 'rb')
    ctg_bases = Counter()
    for i, read in enumerate(bam_fin):
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
    return ctg_covs

  def get_seeds(self):
    if self.options.cheat_seeds:
      self.logger.log('using cheat seeds file: {}'.format(self.options.cheat_seeds))
      seeds = set()
      with open(self.options.cheat_seeds) as fin:
        for line in fin:
          seed = line.strip()
          seeds.add(seed)
      self.logger.log('  - loaded {}'.format(len(seeds)))
      seeds = list(seeds)
    else:
      ctg_size_map = util.get_fasta_sizes(self.options.ctgfasta_path)
      ctg_covs = self.get_ctg_read_counts()
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

  def run(self):
    self.logger.log('get seed contigs from input assembly')

    seeds = [] 
    seeds = self.get_seeds()
    # strip seed contigs into bins such that no more than 4000 bins
    random.shuffle(seeds)
    bins = []
    group_size = max(1, len(seeds) / 4000)
    for i, seed_group in \
      enumerate(util.grouped(seeds, group_size, slop=True)):
      binid = 'bin.{}'.format(i)
      bins.append((binid, seed_group))

    self.logger.log('created {} bins from seeds'.format(len(bins)))
    util.write_pickle(self.options.bins_pickle_path, bins)

