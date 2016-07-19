import os
import pysam
import subprocess
from collections import defaultdict, Counter
import glob
import shutil
from itertools import izip
import random
import numpy as np
from bx.intervals.cluster import ClusterTree

from ..assembler_tools.architect import architect
from ..assembler_tools.barcode_assembler.local import LocalAssembler

from .step import StepChunk
from ..mlib import util
from ..mlib.fq_idx import FastqIndex

MIN_SEED_SIZE = 500
MIN_SEED_SIZE = 400

class AssembleMetaBinnedStep(StepChunk):

  @staticmethod
  def get_steps(options):
    bins = util.load_pickle(options.bins_pickle_path)
    bins2 = util.load_pickle(options.bins2_pickle_path)
    assert len(bins[:-1]) == len(bins2)
    
    for i, ((binid, bcode_set), (_, bcode_counts)) in \
        enumerate(izip(bins[:-1], bins2)):
      if binid != 'bin.contig-100_113':
        continue
      yield AssembleMetaBinnedStep(options, binid)

  def __init__(
    self,
    options,
    binid,
  ):
    self.options = options
    self.binid = binid
    util.mkdir_p(self.outdir)

  def __str__(self):
    return '{}.{}'.format(self.__class__.__name__, self.binid)

  @property
  def outdir(self):
    return self.options.get_bin_dir(self.binid, final=True)

  def outpaths(self, final=False):
    paths = {}
    paths['local-asm-merged.fa'] = os.path.join(self.outdir, 'local-asm-merged.fa')
    #paths['local-asm.p'] = os.path.join(self.outdir, 'local-asm.p')
    #paths['shit'] = 'shit'
    return paths

  def clean_working(self):
    bin_path = self.options.get_bin_dir(self.binid)
    self.logger.log('removing bin directory {}'.format(bin_path))
    shutil.rmtree(bin_path)
    return 

  def run(self):

    self.logger.log('assembling barcoded reads for this bin')
    root_ctg = self.binid[4:]

    ctg_size_map = util.get_fasta_sizes(self.options.ctgfasta_path)

    asmrootdir_path = self.options.get_bin_asm_dir(self.binid)
    util.mkdir_p(asmrootdir_path)

    # create a local assembler for this root contig
    asm = LocalAssembler(
      root_ctg,
      self.options.ctgfasta_path,
      self.options.reads_ctg_bam_path,
      self.options.longranger_fqs_path,
      asmrootdir_path,
      self.logger,
    )

    self.logger.log('determing local assemblies')
    local_asms = asm.gen_local_cands()
    self.logger.log('  - found {} candidates'.format(len(local_asms)))

    # do not locally assemble with other seeds that are lexicographcially
    # smaller than the root. these will get run in other bins
    seed_ctgs = set(filter(
      lambda(c): (
        c > root_ctg and
        ctg_size_map[c] >= MIN_SEED_SIZE
      ),
      ctg_size_map.keys(),
    ))

    self.logger.log('performing local assemblies')
    local_asm_results = asm.assemble(local_asms, filt_ctgs=seed_ctgs)
    self.logger.log('  - finished {}'.format(len(local_asms)))

    #util.write_pickle(
    #  os.path.join(asmrootdir_path, 'local-asm.p'),
    #  asm.dump_info(),
    #)

    # merge output contigs from local assemblies
    self.logger.log('merge long output contigs from local assemblies')
    mergedasm_path = os.path.join(asmrootdir_path, 'local-asm-merged.fa')
    total_asm_contigs = 0
    total_asm_bp = 0
    with open(mergedasm_path, 'w') as fout:
      for i, (local_asm, contig_path) in enumerate(local_asm_results):
        if contig_path == None:
          self.logger.log('contig path for local asm {} not generated'.format(str(local_asm)))
          continue
        fasta = pysam.FastaFile(contig_path)
        for contig in sorted(
          fasta.references,
          key=lambda(c): fasta.get_reference_length(c),
          reverse=True,
        ):
          seq = str(fasta.fetch(contig).upper())
          if len(seq) < 2000:
            break
          total_asm_contigs += 1
          total_asm_bp += len(seq)
          fout.write('>{}.{}${}.{}\n'.format(root_ctg, n_ctg, contig, i))
          fout.write(str(seq) + '\n')

    self.logger.log('  - {} contigs covering {} bases'.format(
      total_asm_contigs,
      total_asm_bp))
      
    def copyfinal(src, dest):
      shutil.copyfile(
        os.path.join(self.options.get_bin_asm_dir(self.binid), src),
        os.path.join(self.options.get_bin_dir(self.binid, final=True), dest),
      )

    self.logger.log('copying deliverables to final')
    copyfinal('local-asm-merged.fa', 'local-asm-merged.fa')
    #copyfinal('local-asm.p', 'local-asm.p')

    self.logger.log('done')

