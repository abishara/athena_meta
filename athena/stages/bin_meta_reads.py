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

  def get_refmapped_seeds(self):
    self.logger.log('determine seeds from ref mappings')

    reffasta_path = '/scratch/PI/serafim/abishara/reference/hmp/all_seqs.fa'
    ctgref_bam_path = '/scratch/users/abishara/scratch/scratch.metagenome/align-contig-hmp.sorted.bam'

    reffasta_path = '/scratch/PI/serafim/abishara/reference/refseq/allref.fa'
    ctgref_bam_path = '/scratch/users/abishara/scratch/scratch.metagenome/align-contig-refseq.sorted.bam'

    # find ref seeds
    self.logger.log('find ref seed contigs')
    ref_seedctg_set = get_ref_seeds(
      reffasta_path,
      ctgref_bam_path,
    )
    # HMP abbrev seeds
    #ref_seedctg_set = set([
    #  'BACT_852|gi|145218786|ref|NZ_AAXE02000112.1|',
    #  # new seeds
    #  #'BACT_847|gi|224462083|gb|ACDQ01000017.1|',
    #  #'BACT_1022|gi|197302419|ref|NZ_ABOU02000032.1|',
    #  #'BACT_847|gi|224462087|gb|ACDQ01000013.1|',
    #  #'BACT_1022|gi|197303624|ref|NZ_ABOU02000050.1|',
    #  'BACT_852|gi|145218676|ref|NZ_AAXE02000002.1|',
    #  # newer seeds
    #  'BACT_545|gi|253795376|ref|NZ_ACOP02000046.1|',
    #  'BACT_496|gi|138274911|ref|NZ_AAXB02000016.1|',
    #  'BACT_852|gi|145218782|ref|NZ_AAXE02000108.1|',
    #  'BACT_171|gi|134301167|ref|NZ_AAVM02000010.1|',
    #])
    refctg_size_map = util.get_fasta_sizes(reffasta_path)
    ref_bases = sum(map(lambda(c): refctg_size_map[c], ref_seedctg_set))
    self.logger.log('  {} seeds covering {} bases'.format(
      len(ref_seedctg_set),
      ref_bases,
    ))

    self.logger.log('find idba0 contigs mapping to ref seeds')
    ref_ctgset_map = get_refmapped_ctgs(
      ref_seedctg_set,
      ctgref_bam_path,
    )
    filtqname_set = set()
    for ref_ctg in ref_ctgset_map:
      filtqname_set |= ref_ctgset_map[ref_ctg]

    fasta = pysam.FastaFile(reffasta_path)
    minirefdir_path = os.path.join(self.options.working_dir, 'miniref')
    util.mkdir_p(minirefdir_path)
    miniref_path = os.path.join(minirefdir_path, 'seedref.fa')
    with open(miniref_path, 'w') as fout:
      for ctg in ref_ctgset_map:
        seq = str(fasta.fetch(ctg).upper())
        fout.write('>{}\n'.format(ctg))
        fout.write('{}\n'.format(seq))
    #util.write_pickle(
    #  os.path.join(self.options.working_dir, 'ref-ctgset.p'),
    #  filtqname_set,
    #)
    #die

    ctg_size_map = util.get_fasta_sizes(self.options.ctgfasta_path)
    self.logger.log('  {} idba0 contigs mapping to seeds'.format(len(filtqname_set)))
    seeds = filter(
      lambda(qname): ctg_size_map[qname] >= MIN_SEED_SIZE,
      filtqname_set,
    )
    self.logger.log('  {} idba0 contigs mapping to seeds  >= {}bp'.format(
      len(seeds),
      MIN_SEED_SIZE,
    ))
    return seeds

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

    #seeds = self.get_refmapped_seeds()
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

#--------------------------------------------------------------------------
# helpers
#--------------------------------------------------------------------------
def get_refmapped_ctgs(
  ref_seed_set,
  ctgref_bam_path,
):

  fin = pysam.Samfile(ctgref_bam_path, 'rb')
  ref_ctgset_map = defaultdict(set)
  for i, read in enumerate(fin):
    if (
      read.is_unmapped or
      read.query_alignment_length < 0.8 * read.query_length
    ):
      continue

    ref_ctg = fin.get_reference_name(read.tid)
    if ref_ctg in ref_seed_set:
      ref_ctgset_map[ref_ctg].add(read.qname)
  fin.close()
  return ref_ctgset_map


def get_ref_seeds(
  reffasta_path,
  ctgref_bam_path,
):

  fin = pysam.Samfile(ctgref_bam_path, 'rb')
  regions_map = defaultdict(lambda: ClusterTree(1,1))
  for i, read in enumerate(fin):
    if read.is_unmapped:
      continue
    ref_ctg = fin.get_reference_name(read.tid)
    regions_map[ref_ctg].insert(read.pos, read.aend, i) 
  fin.close()
  refctg_size_map = util.get_fasta_sizes(reffasta_path)
  fasta = pysam.FastaFile(reffasta_path)
  cand_refctg_set = set()
  cands = []
  for ref_ctg in regions_map.keys():
    bpcov = 0
    for b, e, _ in regions_map[ref_ctg].getregions():
      bpcov += (e - b)

    size = refctg_size_map[ref_ctg]
    cov = 1. * bpcov / size
    cands.append((cov, ref_ctg, size))
    if (
      (size > 50000 and cov > 0.9) or
      (size > 100000 and cov > 0.7)
    ):
      cand_refctg_set.add(ref_ctg)

  for (cov, ref_ctg, size) in sorted(cands, reverse=True)[:20]:
    print 'ref_ctg, cov, size', ref_ctg, cov, size
  return cand_refctg_set

