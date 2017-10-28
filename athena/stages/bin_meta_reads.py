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

random.seed(0)

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
    paths['estimations.p'] = self.options.estimations_path
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
    bam_fin.close()
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

    # estimate fragment lengths from mapping to seeds
    estimated_len = self.get_estimated_frag_length()
    util.write_pickle(self.options.estimations_path, estimated_len)

  def get_estimated_frag_length(self):

    self.logger.log('estimating extracted DNA fragment size from alignments to large contigs')
    MIN_CTG_SIZE = 50000
    ctg_size_map = util.get_fasta_sizes(self.options.ctgfasta_path)
    # estimate if at least 1mb of sequence in contigs of threshold size
    ctg_sizes = filter(lambda(c, l): l >= MIN_CTG_SIZE, ctg_size_map.items())
    substrate_size = sum(map(lambda(c,v): v, ctg_sizes))
    if substrate_size < 1e06:
      self.logger.log('  - <1Mb of sequence in contigs >={}bp'.format(MIN_CTG_SIZE))
      self.logger.log('  - defaulting to 10kb estimate')
      return 10000

    # randomly subsample sufficiently large contigs to prevent from using all
    # the DNA from one or two organisms to be used for estimation
    random.shuffle(ctg_sizes)
    substrate_size = 0
    sample_ctgs = set()
    for ctg, size in ctg_sizes:
      substrate_size += size
      if substrate_size > 1e06:
        break
      sample_ctgs.add(ctg)

    fragment_sizes = []
    bam_fin = pysam.Samfile(self.options.reads_ctg_bam_path, 'rb')
    for ctg in sample_ctgs:
      for (bcode, reads_iter) in get_bcode_iter(bam_fin.fetch(ctg)):
        if bcode == None:
          continue
        for reads in get_cloud_iter(reads_iter, bam_fin):
          num_reads = len(reads)
          fragment_size = (
            max(reads,key=lambda(r): r.aend).aend -
            min(reads,key=lambda(r): r.pos).pos
          )
          if num_reads >= 6:
            fragment_sizes.append(fragment_size)
    bam_fin.close()
    estimated_len = int(np.median(fragment_sizes))
    self.logger.log('  - median fragment length of {}bp from {} imputed long fragments mapped to long seed contigs'.format(
      estimated_len, len(fragment_sizes)))
    return estimated_len
      
def get_bcode_iter(bam_fin):
  seen_set = set()
  for bcode, reads_iter in groupby(
    sorted(
      bam_fin,
      key=lambda(read): (util.get_barcode(read), read.pos),
    ),
    key=lambda(read): util.get_barcode(read),
  ):
    if bcode != None:
      yield bcode, reads_iter
    assert bcode not in seen_set, \
      "input bam {} not in bcode sorted order".format(bam_path)
  raise StopIteration

def get_cloud_iter(reads_iter, bam_fin):
  MIN_CLOUD_GAP = 20000
  unmapped_reads = set()
  skey_func = lambda(r): (r.is_unmapped, r.tid, r.pos)
  gkey_func = lambda(r): (r.is_unmapped, r.tid)
  for (umap, tid), treads_iter in groupby(
    sorted(reads_iter, key=skey_func), key=gkey_func,
  ):
    if umap:
      continue
    creads = []
    laend = None
    for read in treads_iter:
      if laend == None or read.pos - laend < MIN_CLOUD_GAP:
        creads.append(read)
        laend = read.aend
      else:
        assert read.pos >= laend + MIN_CLOUD_GAP
        yield creads
        creads = [read]
        laend = None
    if len(creads) > 1:
      yield creads
  raise StopIteration

