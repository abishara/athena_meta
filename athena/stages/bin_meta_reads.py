import os
import pysam
import shutil
from bx.intervals.cluster import ClusterTree
from collections import defaultdict, Counter
import numpy as np
import networkx as nx
from itertools import groupby, combinations
import random

from .step import StepChunk
from ..mlib import util

MIN_SEED_SIZE = 500
MIN_SEED_SIZE = 400

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

  def get_seeds(self):
    ctg_size_map = util.get_fasta_sizes(self.options.ctgfasta_path)
    seeds = ctg_size_map.keys()
    seeds = filter(
      lambda(qname): ctg_size_map[qname] >= MIN_SEED_SIZE,
      seeds,
    )
    self.logger.log('  {} idba0 contigs >= {}bp'.format(
      len(seeds),
      MIN_SEED_SIZE,
    ))
    return seeds

  def run(self):
    self.logger.log('get seed contigs from idba0')

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

