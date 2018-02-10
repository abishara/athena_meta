import os
import pysam
import subprocess
import shutil
from collections import defaultdict, Counter

from .step import StepChunk
from ..mlib import util

# NOTE must be in path
flyebin_path = 'flye'

class AssembleOLCStep(StepChunk):

  @staticmethod
  def get_steps(options):
    assert os.path.isfile(options.ctgfasta_path+'.sa'), \
      "ctgfasta_path option in config must reference a BWA index"
    yield AssembleOLCStep(options)

  @staticmethod
  def deliver_message(options):
    return "Athena contigs: {}".format(os.path.join(
      options.results_dir,
      'olc',
      'athena.asm.fa',
    ))

  @property
  def outdir(self):
    return os.path.join(self.options.results_dir, 'olc')

  def outpaths(self, final=False):
    paths = {}
    paths['athena.asm.fa'] = os.path.join(self.outdir, 'athena.asm.fa')
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
    self.logger.log('jointly overlap-assemble')


    premergedfa_path = os.path.join(self.outdir, 'pre-flye-input-contigs.fa')
    premergedfiltfa_path = os.path.join(self.outdir, 'pre-flye-input-contigs.filt.fa')
    seedsfa_path = os.path.join(self.outdir, 'seed-contigs.fa')
    mergedfiltfa_path = os.path.join(self.outdir, 'flye-input-contigs.fa')

    # load all the bins and merge
    if not os.path.isfile(premergedfa_path):
      self.logger.log('merge input contigs')
      bins = util.load_pickle(self.options.bins_pickle_path)
      input_paths = []
      seed_ctgs = set()
      for (binid, seed_group) in bins:
        for ctg in seed_group:
          seed_ctgs.add(ctg)
        bindir_path = self.options.get_bin_dir(binid, final=True)
        fa_path = os.path.join(bindir_path, 'local-asm-merged.fa')
        if not os.path.isfile(fa_path):
          self.logger.log(' not found, skipping {}'.format(fa_path))
        elif not is_valid_fasta(fa_path):
          self.logger.log('WARNING input fasta not valid: {}'.format(fa_path))
        else:
          input_paths.append(fa_path)
      util.concat_files(input_paths, premergedfa_path)
    assert is_valid_fasta(premergedfa_path), "merge FASTA not valid"

    # filter subassembled inputs that are not useful
    mergedbam_path = os.path.join(self.outdir, 'align-inputs.bam')
    cmd = 'bwa mem -t 4 {} {} | samtools view -bS - | samtools sort -o {} - '.format(
      self.options.ctgfasta_path,
      premergedfa_path,
      mergedbam_path,
    )
    if not os.path.isfile(mergedbam_path):
      print 'cmd', cmd
      subprocess.check_call(cmd, shell=True)
      cmd = 'samtools index {}'.format(mergedbam_path)
      print 'cmd', cmd
      subprocess.check_call(cmd, shell=True)

    if not os.path.isfile(mergedfiltfa_path):
      self.logger.log('filter short subassembled contigs and merge with seeds')
      filter_inputs(
        mergedbam_path,
        premergedfa_path,
        premergedfiltfa_path,
      )

      # append input seed contigs
      with open(seedsfa_path, 'w') as fout:
        fasta = pysam.FastaFile(self.options.ctgfasta_path)
        seed_ctgs = self.load_seeds()
        for ctg in seed_ctgs:  
          size = fasta.get_reference_length(ctg)
          if size < 1000:
            continue
          seq = str(fasta.fetch(ctg).upper())
          for i in xrange(5):
          #for i in xrange(2):
            fout.write('>{}.{}\n'.format(ctg, i))
            fout.write(str(seq) + '\n')

      util.concat_files(
        [premergedfiltfa_path, seedsfa_path], 
        mergedfiltfa_path,
      )
    assert is_valid_fasta(mergedfiltfa_path), "merge FASTA not valid"

    # run flye OLC assembly
    seed_draft_size = 0
    fasta = pysam.FastaFile(self.options.ctgfasta_path)
    seed_ctgs = self.load_seeds()
    for ctg in seed_ctgs:
      seed_draft_size += fasta.get_reference_length(ctg)

    flye0_path = os.path.join(self.outdir, 'flye-asm-1')
    flye_contigs_path = os.path.join(flye0_path, 'scaffolds.fasta')
    cmd = '{} --subassemblies {} --out-dir {} --genome-size {} --threads 4 --min-overlap 1000'.format(
      flyebin_path,
      mergedfiltfa_path,
      flye0_path,
      seed_draft_size,
    )
    if not os.path.isfile(flye_contigs_path):
      print 'launching Flye OLC assembly'
      print 'cmd', cmd
      subprocess.check_call(cmd, shell=True)

    # copy flye output as athena final output
    final_fa_path = os.path.join(self.outdir, 'athena.asm.fa')
    shutil.copy(flye_contigs_path, final_fa_path)
    self.logger.log('done')

#--------------------------------------------------------------------------
# helpers
#--------------------------------------------------------------------------
  def load_seeds(self):
    bins = util.load_pickle(self.options.bins_pickle_path)
    seed_ctgs = set()
    for (binid, seed_group) in bins:
      for ctg in seed_group:
        seed_ctgs.add(ctg)
    return seed_ctgs

def filter_inputs(
  mergedbam_path,
  mergedfa_path,
  mergedfiltfa_path,
):
  ctg_size_map = util.get_fasta_sizes(mergedfa_path)
  fhandle = pysam.Samfile(mergedbam_path)
  full_ctgs = set()
  for read in fhandle:
    if read.is_unmapped:
      continue
    if read.query_alignment_length + 1000 >= ctg_size_map[read.qname]:
      full_ctgs.add(read.qname)
  fhandle.close()

  print 'orig ctgs', len(ctg_size_map)
  new_ctgs = set(ctg_size_map.keys()) - full_ctgs
  print 'filtered ctgs', len(new_ctgs)
  fasta = pysam.FastaFile(mergedfa_path)
  with open(mergedfiltfa_path, 'w') as fout:
    for ctg in new_ctgs:
      seq = str(fasta.fetch(ctg).upper())
      fout.write('>{}\n'.format(ctg))
      fout.write('{}\n'.format(seq))

def is_valid_fasta(fa_path):
  alpha = set('acgtACGT')
  with open(fa_path) as fin:
    for line in fin:
      if line.startswith('>'):
        continue
      seq = set(line.strip())
      if not seq.issubset(alpha):
        return False

  return True

