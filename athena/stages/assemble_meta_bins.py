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
    
    for i, binid in bins:
      #if binid != 'bin.contig-100_113':
      #  continue
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

    def get_barcode(read):
      filt_list = filter(lambda(k, v): k == 'BX', read.tags)
      if filt_list == []: 
        return None
      else:
        k, v = filt_list[0]
        return v

    # check enough read/barcode coverage to warrant denovo assembly
    numreads = 0
    bcodes = set()
    fhandle = pysam.Samfile(self.options.reads_ctg_bam_path, 'rb')
    for read in fhandle.fetch(root_ctg):
      if read.is_unmapped or read.mapq < 10:
        continue
      numreads += 1
      bcode = get_barcode(read)
      if bcode != None:
        bcodes.add(bcode)
    fhandle.close()
    size = ctg_size_map[ctg]
    cov = 95. * numreads / size
    if cov < 10. or len(bcodes) < 30.:
      self.logger.log('seed contig does not have high enough coverage')
      self.logger.log('  - {} bcodes, {}x'.format(len(bcodes), cov))
      util.touch(self.outpaths()['local-asm-merged.fa'])
      self.logger.log('done')
      return

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
          fout.write('>{}.{}${}.{}\n'.format(local_asm.root_ctg, local_asm.link_ctg, contig, i))
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

