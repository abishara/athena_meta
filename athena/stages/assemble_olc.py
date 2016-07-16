import os
import pysam
import subprocess
import shutil
from collections import defaultdict, Counter

from .step import StepChunk
from ..mlib import util

# FIXME harcoded paths
# NOTE must be in path
canubin_path = '/home/abishara/sources/canu/Linux-amd64/bin/canu'
olapbin_path = '/home/abishara/sources/canu/Linux-amd64/bin/ovStoreDump'

class AssembleOLCStep(StepChunk):

  @staticmethod
  def get_steps(options):
    yield AssembleOLCStep(options)

  @property
  def outdir(self):
    return os.path.join(self.options.results_dir, 'olc')

  def outpaths(self, final=False):
    paths = {}
    paths['shit'] = 'shit'
    return paths

  def __init__(
    self,
    options,
  ):
    self.options = options
    util.mkdir_p(self.outdir)

  def __str__(self):
    return self.__class__.__name__

  def run(self):
    self.logger.log('jointly assemble bins with OLC and architect scaffolding')

    # collect input contig and local reasm contigs
    self.logger.log('merge input contigs')

    mergedfa_path = os.path.join(self.outdir, 'canu-input-contigs.fa')
    seedsfa_path = os.path.join(self.outdir, 'hmp-seed-contigs.fa')

    # load all the bins
    bins = util.load_pickle(self.options.bins_pickle_path)
    # skip union bin
    bins = bins[:-1]

    #input_paths = [self.options.ctgfasta_path]
    input_paths = []
    seed_ctgs = set()
    for (binid, _) in  bins:
      seed_ctg = binid[4:]
      seed_ctgs.add(seed_ctg)
      bindir_path = self.options.get_bin_dir(binid, final=True)
      fa_path = os.path.join(bindir_path, 'local-asm-merged.fa')
      if not os.path.isfile(fa_path):
        self.logger.log(' not found, skipping {}'.format(fa_path))
      else:
        input_paths.append(fa_path)

    # append input hmp seed contigs
    with open(seedsfa_path, 'w') as fout:
      fasta = pysam.FastaFile(self.options.ctgfasta_path)
      for ctg in seed_ctgs:  
        seq = str(fasta.fetch(ctg).upper())
        fout.write('>{}\n'.format(ctg))
        fout.write(str(seq) + '\n')
    input_paths.append(seedsfa_path)

    #util.concat_files(input_paths, mergedfa_path)
    #die

    #self.logger.log('  {} contigs, covering {} bases'.format(
    #  total_asm_contigs,
    #  total_asm_bp,
    #))
    with util.cd(self.outdir):
      canu0_path = 'canu-asm-0'
      mergedfa_path = 'canu-input-contigs.fa'

      cmd = '{} useGrid=0 errorRate=0.07 genomeSize=5.00m stopOnReadQuality=false -d {} -p canu -pacbio-corrected {}'.format(
        canubin_path,
        canu0_path,
        mergedfa_path
      )
      print 'cmd', cmd
      subprocess.check_call(cmd, shell=True)

    self.logger.log('done')

#--------------------------------------------------------------------------
# helpers
#--------------------------------------------------------------------------

