import os
import pysam
import shutil
from bx.intervals.cluster import ClusterTree
from collections import defaultdict, Counter
import numpy as np
import networkx as nx
from itertools import groupby, combinations
import community

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
    paths['bins2.p'] = self.options.bins2_pickle_path
    #paths['shit.p'] = 'shit'
    return paths

  def run(self):
    self.logger.log('determine bins from seed contigs')

    ctg_size_map = util.get_fasta_sizes(self.options.ctgfasta_path)
    
    hmpfasta_path = '/scratch/PI/serafim/abishara/reference/hmp/all_seqs.fa'
    ctghmp_bam_path = '/scratch/users/abishara/scratch/scratch.metagenome/align-contig-hmp.sorted.bam'

    # find hmp seeds
    self.logger.log('find hmp seed contigs')
    #hmp_seedctg_set = get_hmp_seeds(
    #  hmpfasta_path,
    #  ctghmp_bam_path,
    #)
    hmp_seedctg_set = set([
      'BACT_852|gi|145218786|ref|NZ_AAXE02000112.1|',
      # new seeds
      #'BACT_847|gi|224462083|gb|ACDQ01000017.1|',
      #'BACT_1022|gi|197302419|ref|NZ_ABOU02000032.1|',
      #'BACT_847|gi|224462087|gb|ACDQ01000013.1|',
      #'BACT_1022|gi|197303624|ref|NZ_ABOU02000050.1|',
      'BACT_852|gi|145218676|ref|NZ_AAXE02000002.1|',
      # newer seeds
      'BACT_545|gi|253795376|ref|NZ_ACOP02000046.1|',
      'BACT_496|gi|138274911|ref|NZ_AAXB02000016.1|',
      'BACT_852|gi|145218782|ref|NZ_AAXE02000108.1|',
      'BACT_171|gi|134301167|ref|NZ_AAVM02000010.1|',
    ])
    hmpctg_size_map = util.get_fasta_sizes(hmpfasta_path)
    hmp_bases = sum(map(lambda(c): hmpctg_size_map[c], hmp_seedctg_set))
    self.logger.log('  {} seeds covering {} bases'.format(
      len(hmp_seedctg_set),
      hmp_bases,
    ))

    self.logger.log('find idba0 contigs mapping to hmp seeds')
    hmp_ctgset_map = get_hmpmapped_ctgs(
      hmp_seedctg_set,
      ctghmp_bam_path,
    )
    filtqname_set = set()
    for hmp_ctg in hmp_ctgset_map:
      filtqname_set |= hmp_ctgset_map[hmp_ctg]

    fasta = pysam.FastaFile(hmpfasta_path)
    minirefdir_path = os.path.join(self.options.working_dir, 'miniref')
    util.mkdir_p(minirefdir_path)
    miniref_path = os.path.join(minirefdir_path, 'hmpseedref.fa')
    with open(miniref_path, 'w') as fout:
      for ctg in hmp_ctgset_map:
        seq = str(fasta.fetch(ctg).upper())
        fout.write('>{}\n'.format(ctg))
        fout.write('{}\n'.format(seq))
    #util.write_pickle(
    #  os.path.join(self.options.working_dir, 'hmp-ctgset.p'),
    #  filtqname_set,
    #)
    #die

    self.logger.log('  {} idba0 contigs mapping to seeds'.format(len(filtqname_set)))

    tmp_hits_path = '/home/abishara/tmp/tmp_bcode-ctg-hits.txt'
    hits_path      = self.options.bcode_ctg_hits_path + '_bcode-ctg-hits.txt'

    bcode_idx_map = self.options.bcode_idx_map
    idx_bcode_map = {v: k for k, v in bcode_idx_map.items()}
    ctg_idx_map = self.options.ctg_idx_map
    idx_ctg_map = {v: k for k, v in ctg_idx_map.items()}

    # some contigs may have no barcode hits...
    filtqname_set = set(filter(
      lambda(qname): qname in ctg_idx_map,
      filtqname_set,
    ))
    self.logger.log('  {} idba0 contigs mapping to seeds with barcode hits'.format(len(filtqname_set)))
    filtqname_set = set(filter(
      lambda(qname): ctg_size_map[qname] >= MIN_SEED_SIZE,
      filtqname_set,
    ))
    self.logger.log('  {} idba0 contigs mapping to seeds with barcode hits >= {}bp'.format(
      len(filtqname_set),
      MIN_SEED_SIZE,
    ))
    filt_ctgidx_set = set(map(
      lambda(qname): ctg_idx_map[qname],
      filtqname_set,
    ))

    self.logger.log('parsing barcode-contig hits to determine barcode sets to compile')
    ctg_bcodes_map = defaultdict(Counter)
    total_reads = 0
    #fout = open(tmp_hits_path, 'w')
    with open(tmp_hits_path) as f:
    #with open(hits_path) as f:
      for k, lines in groupby(f, lambda(l): l.split('\t')[0]):
        hits_set = set()
        for line in lines:
          bcode_idx, ctg_idx, cnt = map(int, line.split('\t'))
          # skip irrelevant contig hits
          if ctg_idx not in filt_ctgidx_set:
            continue
          total_reads += cnt
          hits_set.add(ctg_idx)
          ctg_bcodes_map[ctg_idx][bcode_idx] += cnt
          #fout.write(line)
    #fout.close()
    
    #print 'total reads', total_reads

    # create a bin for each seed contig with sufficient coverage
    bins = []
    bins2 = []
    skipped = 0
    all_bcode_set = set()
    for i, ctg in enumerate(filtqname_set):
      ctg_idx = ctg_idx_map[ctg]
      bcode_counts = ctg_bcodes_map[ctg_idx]
      #bcode_set = set(bcode_counts.keys())
      bcode_set = set(map(lambda(x): idx_bcode_map[x], bcode_counts.keys()))
      numreads = sum(bcode_counts.values())
      size = ctg_size_map[ctg]
      cov = 95. * numreads / size
      if cov < 10. or len(bcode_set) < 30.:
      #if cov < 20. or len(bcode_set) < 100:
        self.logger.log(' skipping {} for cov: {} depth, {} bcodes'.format(
          ctg,
          cov,
          len(bcode_set)
        ))
        skipped += 1
        continue
      all_bcode_set |= bcode_set
      binid = 'bin.{}'.format(ctg)

      bins.append((binid, bcode_set))
      bins2.append((binid, bcode_counts))
    bins.append(('bin.union', all_bcode_set))
      
    self.logger.log('created {} bins from seeds'.format(len(bins)))
    self.logger.log('  - {} skipped for coverage'.format(skipped))
    util.write_pickle(self.options.bins_pickle_path, bins)
    util.write_pickle(self.options.bins2_pickle_path, bins2)

#--------------------------------------------------------------------------
# helpers
#--------------------------------------------------------------------------
def get_hmpmapped_ctgs(
  hmp_seed_set,
  ctghmp_bam_path,
):

  fin = pysam.Samfile(ctghmp_bam_path, 'rb')
  hmp_ctgset_map = defaultdict(set)
  for i, read in enumerate(fin):
    if (
      read.is_unmapped or
      read.query_alignment_length < 0.8 * read.query_length
    ):
      continue

    hmp_ctg = fin.get_reference_name(read.tid)
    if hmp_ctg in hmp_seed_set:
      hmp_ctgset_map[hmp_ctg].add(read.qname)
  fin.close()
  return hmp_ctgset_map


def get_hmp_seeds(
  hmpfasta_path,
  ctghmp_bam_path,
):

  fin = pysam.Samfile(ctghmp_bam_path, 'rb')
  regions_map = defaultdict(lambda: ClusterTree(1,1))
  for i, read in enumerate(fin):
    if read.is_unmapped:
      continue
    hmp_ctg = fin.get_reference_name(read.tid)
    regions_map[hmp_ctg].insert(read.pos, read.aend, i)
  fin.close()

  hmpctg_size_map = util.get_fasta_sizes(hmpfasta_path)
  fasta = pysam.FastaFile(hmpfasta_path)
  cand_hmpctg_set = set()
  for hmp_ctg in regions_map.keys():
    bpcov = 0
    for b, e, _ in regions_map[hmp_ctg].getregions():
      bpcov += (e - b)

    size = hmpctg_size_map[hmp_ctg]
    cov = 1. * bpcov / size
    if size > 50000 and cov > 0.9:
      cand_hmpctg_set.add(hmp_ctg)

  return cand_hmpctg_set

