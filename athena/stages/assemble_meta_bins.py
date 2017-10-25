import os
import pysam
import subprocess
from collections import defaultdict, Counter
import shutil

from ..assembler_tools.barcode_assembler.local import LocalAssembler
from ..assembler_tools.barcode_assembler import local as local_assembler

from .step import StepChunk
from ..mlib import util
from ..mlib.fq_idx import FastqIndex

MIN_SEED_SIZE = 500
MIN_SEED_SIZE = 400

class AssembleMetaBinnedStep(StepChunk):

  @staticmethod
  def get_steps(options):
    util.check_prereqs()
    bins = util.load_pickle(options.bins_pickle_path)
    
    for i, (binid, seeds) in enumerate(bins):
      yield AssembleMetaBinnedStep(options, binid, seeds)

  def __init__(
    self,
    options,
    binid,
    seeds,
  ):
    self.options = options
    self.binid = binid
    self.seeds = seeds
    util.mkdir_p(self.outdir)

  def __str__(self):
    return '{}.{}'.format(self.__class__.__name__, self.binid)

  @property
  def outdir(self):
    return self.options.get_bin_dir(self.binid, final=True)

  def outpaths(self, final=False):
    paths = {}
    paths['local-asm-merged.fa'] = os.path.join(self.outdir, 'local-asm-merged.fa')
    #paths['shit'] = 'shit'
    return paths

  def clean_working(self):
    bin_path = self.options.get_bin_dir(self.binid)
    self.logger.log('removing bin directory {}'.format(bin_path))
    shutil.rmtree(bin_path)
    return 

  def run(self):

    self.logger.log('performing local assembly for {} seeds'.format(
      len(self.seeds)))

    local_assembler.DS_SUBASM_COV = self.options.ds_subasm_cov
    local_assembler.SEED_SELF_ASM_SIZE = self.options.seed_self_asm_size
    self.logger.log('targeting {}x short-read subassembly coverage'.format(
      local_assembler.DS_SUBASM_COV,
    ))
    self.logger.log('using barcodes mapped within {}bp from seed end-points for seed subassembly'.format(
      local_assembler.SEED_SELF_ASM_SIZE,
    ))

    asmrootdir_path = self.options.get_bin_asm_dir(self.binid)
    util.mkdir_p(asmrootdir_path)

    out_paths = []
    for ctg in self.seeds:
      asmdir = os.path.join(asmrootdir_path, ctg)
      util.mkdir_p(asmdir)
      pass_path = os.path.join(asmdir, 'pass')
      if os.path.isfile(pass_path):
        self.logger.log('seed {} output already generated'.format(ctg))
        continue
      self.do_local_assembly(ctg, asmdir)
      out_path = os.path.join(asmdir, 'local-asm-merged.fa')
      assert os.path.isfile(out_path), \
        "output does not exist {}".format(out_path)
      out_paths.append(out_path)

    # concatenate all results from seed contigs
    self.logger.log('merging all outputs')
    mergedout_path = os.path.join(asmrootdir_path, 'local-asm-merged.fa')
    util.concat_files(out_paths, mergedout_path)
    shutil.copy(mergedout_path, self.options.get_bin_dir(self.binid, final=True))
    self.logger.log('done')

  def do_local_assembly(self, root_ctg, asmrootdir_path):

    self.logger.log('assembling barcoded reads for seed {}'.format(root_ctg))

    ctg_size_map = util.get_fasta_sizes(self.options.ctgfasta_path)

    # check enough read/barcode coverage to warrant denovo assembly
    numreads = 0
    bcodes = set()
    fhandle = pysam.Samfile(self.options.reads_ctg_bam_path, 'rb')
    for read in fhandle.fetch(root_ctg):
      if read.is_unmapped or read.mapq < 10:
        continue
      numreads += 1
      bcode = util.get_barcode(read)
      if bcode != None:
        bcodes.add(bcode)
    fhandle.close()
    size = ctg_size_map[root_ctg]
    cov = 95. * numreads / size
    if cov < 10. or len(bcodes) < 30.:
      self.logger.log('seed {} contig does not have high enough coverage'.format(root_ctg))
      self.logger.log('  - {} bcodes, {}x'.format(len(bcodes), cov))
      out_path = os.path.join(asmrootdir_path, 'local-asm-merged.fa')
      util.touch(out_path)
      return

    # create a local assembler for this root contig
    asm = LocalAssembler(
      root_ctg,
      self.options.ctgfasta_path,
      self.options.reads_ctg_bam_path,
      self.options.input_fqs,
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
          link_name = local_asm.link_ctg
          if link_name == None:
            link_name = 'seed'
          fout.write('>{}.{}${}.{}\n'.format(local_asm.root_ctg, link_name, contig, i))
          fout.write(str(seq) + '\n')

    self.logger.log('  - {} contigs covering {} bases'.format(
      total_asm_contigs,
      total_asm_bp))
    pass_path = os.path.join(asmrootdir_path, 'pass')
    util.touch(pass_path)
      
