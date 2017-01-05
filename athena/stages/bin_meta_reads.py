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

    #reffasta_path = '/scratch/PI/serafim/abishara/reference/refseq/allref.fa'
    #ctgref_bam_path = '/scratch/users/abishara/scratch/scratch.metagenome/align-contig-refseq.sorted.bam'

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
    #seeds = self.get_seeds()
    seeds = [
'NODE_24_length_73415_cov_20.5401_ID_211432',
'NODE_1012_length_1696_cov_4.24924_ID_213452',
'NODE_20_length_74840_cov_18.113_ID_211424',
'NODE_899_length_1972_cov_4.4554_ID_213222',
'NODE_195_length_8296_cov_18.4116_ID_211806',
'NODE_26_length_65931_cov_20.7105_ID_211444',
'NODE_9218_length_225_cov_0.905882_ID_229822',
'NODE_11814_length_144_cov_1.96629_ID_235000',
'NODE_11691_length_189_cov_12.3582_ID_234754',
'NODE_42_length_52586_cov_19.0512_ID_211470',
'NODE_7_length_148007_cov_24.7768_ID_211386',
'NODE_10561_length_214_cov_1.01258_ID_232500',
'NODE_715_length_2547_cov_3.34952_ID_212850',
'NODE_12332_length_107_cov_18.5769_ID_236036',
'NODE_1226_length_1336_cov_317.912_ID_213890',
'NODE_45_length_47671_cov_19.7511_ID_211480',
'NODE_23_length_73897_cov_21.4318_ID_211430',
'NODE_13575_length_86_cov_9.41935_ID_238522',
'NODE_13575_length_86_cov_9.41935_ID_238522',
'NODE_8181_length_235_cov_1.07778_ID_227760',
'NODE_11348_length_208_cov_1_ID_234070',
'NODE_115_length_16070_cov_28.8975_ID_211636',
'NODE_12131_length_113_cov_11.1034_ID_235634',
'NODE_28_length_64627_cov_23.0826_ID_211448',
'NODE_6_length_172776_cov_23.2538_ID_211384',
'NODE_9602_length_222_cov_1.57485_ID_230586',
'NODE_112_length_16940_cov_17.0784_ID_211630',
'NODE_5658_length_265_cov_1.1381_ID_222738',
'NODE_841_length_2135_cov_51.8745_ID_213108',
'NODE_191_length_8424_cov_21.3788_ID_211800',
'NODE_12332_length_107_cov_18.5769_ID_236036',
'NODE_201_length_7899_cov_18.6523_ID_211822',
'NODE_99_length_21438_cov_26.2619_ID_211604',
'NODE_1354_length_1142_cov_19.8096_ID_214142',
'NODE_111_length_17433_cov_14.0972_ID_211628',
'NODE_298_length_5672_cov_4.05875_ID_212006',
'NODE_1088_length_1549_cov_26.6178_ID_213614',
'NODE_3200_length_381_cov_13.862_ID_217872',
'NODE_3201_length_381_cov_19.2822_ID_217874',
'NODE_131_length_14453_cov_24.2464_ID_211668',
'NODE_21_length_74503_cov_24.8564_ID_211426',
'NODE_21_length_74503_cov_24.8564_ID_211426',
'NODE_106_length_19326_cov_19.4395_ID_211618',
'NODE_9_length_119506_cov_24.3114_ID_211390',
'NODE_6996_length_248_cov_1.51295_ID_225400',
'NODE_6375_length_255_cov_1.225_ID_224168',
'NODE_10944_length_211_cov_1_ID_233262',
'NODE_2_length_240395_cov_25.0891_ID_211376',
'NODE_11808_length_146_cov_1.9011_ID_234988',
'NODE_10589_length_214_cov_1.02516_ID_232556',
'NODE_11924_length_131_cov_84.9342_ID_235220',
'NODE_10472_length_215_cov_1.1125_ID_232322',
'NODE_2405_length_512_cov_26.0197_ID_216278',
'NODE_1_length_290890_cov_20.7303_ID_211374',
'NODE_7536_length_242_cov_1.27807_ID_226476',
'NODE_11707_length_174_cov_1.36134_ID_234786',
'NODE_9597_length_222_cov_1.07186_ID_230576',
'NODE_7345_length_244_cov_1.07407_ID_226098',
'NODE_6454_length_254_cov_1.35176_ID_224326',
'NODE_9205_length_225_cov_1.00588_ID_229796',
    ]

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
      (size > 100000 and cov > 0.9) or
      (size > 300000 and cov > 0.7)
    ):
    #if (
    #  (size > 10000 and cov > 0.9) or
    #  (size > 100000 and cov > 0.7)
    #):
      cand_refctg_set.add(ref_ctg)

  for (cov, ref_ctg, size) in sorted(cands, reverse=True)[:20]:
    print 'ref_ctg, cov, size', ref_ctg, cov, size
  return cand_refctg_set

