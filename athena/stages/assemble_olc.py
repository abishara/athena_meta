import os
import pysam
import subprocess
import shutil
from collections import defaultdict, Counter

from .step import StepChunk
from ..mlib import util

# NOTE must be in path
canubin_path = 'canu'

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
    self.logger.log('jointly assemble bins with OLC')

    # collect input contig and local reasm contigs
    self.logger.log('merge input contigs')

    premergedfa_path = os.path.join(self.outdir, 'pre-canu-input-contigs.fa')
    premergedfiltfa_path = os.path.join(self.outdir, 'pre-canu-input-contigs.filt.fa')
    seedsfa_path = os.path.join(self.outdir, 'seed-contigs.fa')
    mergedfiltfa_path = os.path.join(self.outdir, 'canu-input-contigs.fa')

    # load all the bins
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

    if not os.path.isfile(premergedfa_path):
      util.concat_files(input_paths, premergedfa_path)
    assert is_valid_fasta(premergedfa_path), "merge FASTA not valid"

    mergedbam_path = os.path.join(self.outdir, 'align-inputs.bam')
    cmd = 'bwa mem -t 8 {} {} | samtools view -bS - | samtools sort -o {} - '.format(
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

    filter_inputs(
      mergedbam_path,
      premergedfa_path,
      premergedfiltfa_path,
    )

    # append input seed contigs
    with open(seedsfa_path, 'w') as fout:
      fasta = pysam.FastaFile(self.options.ctgfasta_path)
      for ctg in seed_ctgs:  
        seq = str(fasta.fetch(ctg).upper())
        for i in xrange(5):
        #for i in xrange(2):
          fout.write('>{}.{}\n'.format(ctg, i))
          fout.write(str(seq) + '\n')

    if not os.path.isfile(mergedfiltfa_path):
      util.concat_files(
        [premergedfiltfa_path, seedsfa_path], 
        mergedfiltfa_path,
      )
    assert is_valid_fasta(mergedfiltfa_path), "merge FASTA not valid"

    canu0_path = os.path.join(self.outdir, 'canu-asm-1')
    cmd = \
'{} \
gnuplotTested=true \
useGrid=0  \
correctedErrorRate=0.06  \
genomeSize=45.00m  \
contigFilter="2 2000 1.0 1.0 2" \
stopOnReadQuality=false \
-d {} \
-p canu \
-pacbio-corrected {}'.format(
      canubin_path,
      canu0_path,
      mergedfiltfa_path
    )
# stale unused canu flags
#gridOptions="-p owners" \
#oeaMemory=12 cnsMemory=32 batMemory=50 \
     #die
    canu_contigs_path = os.path.join(canu0_path, 'canu.contigs.fasta')
    if not os.path.isfile(canu_contigs_path):
      print 'launching OLC assembly'
      print 'cmd', cmd
      subprocess.check_call(cmd, shell=True)

    # index assembled contigs 
    self.logger.log('index canu assembled contigs')
    cmd = 'bwa index {}'.format(canu_contigs_path)
    subprocess.check_call(cmd, shell=True)

    # align init contigs to canu contigs
    self.logger.log('aligning seed contigs to olc contigs')
    inputfa_path = self.options.ctgfasta_path
    aligninputsbam_path = os.path.join(self.outdir, 'align-input-contigs.canu-contigs.bam')
    cmd = 'bwa mem -t 8 {} {} | samtools view -bS - | samtools sort -o {} -'.format(
      canu_contigs_path,
      inputfa_path,
      aligninputsbam_path,
    )
    print 'cmd', cmd
    subprocess.check_call(cmd, shell=True)
    cmd = 'samtools index {}'.format(aligninputsbam_path)
    print 'cmd', cmd
    subprocess.check_call(cmd, shell=True)

    seeds = set() 
    bins = util.load_pickle(self.options.bins_pickle_path)
    for _, _seeds in bins:
      seeds |= set(_seeds)

    input_fasta = pysam.FastaFile(inputfa_path)
    self.logger.log('determine unmapped contigs out of {} input contigs'.format(
      input_fasta.nreferences))

    unmap_input_ctgs = get_unmapped_ctgs(
      inputfa_path,
      aligninputsbam_path,
    )
    self.logger.log('  - {} unmapped'.format(len(unmap_input_ctgs)))
    self.logger.log('  - {}/{} seeds unmapped'.format(
      len(unmap_input_ctgs & seeds),
      len(seeds),
    ))
    input_unmapfa_path = os.path.join(self.outdir, 'input.unmap.fa')
    with open(input_unmapfa_path, 'w') as fout:
      for ctg in unmap_input_ctgs:
        seq = str(input_fasta.fetch(ctg).upper())
        fout.write('>{}\n'.format(ctg))
        fout.write('{}\n'.format(seq))

    final_fa_path = os.path.join(self.outdir, 'athena.asm.fa')
    util.concat_files(
      [input_unmapfa_path, canu_contigs_path],
      final_fa_path,
    )

    self.logger.log('done')

#--------------------------------------------------------------------------
# helpers
#--------------------------------------------------------------------------
def get_unmapped_ctgs(fa_path, bam_path):

  ctg_size_map = util.get_fasta_sizes(fa_path)
  fhandle = pysam.Samfile(bam_path, 'rb')
  map_ctgs = set()
  for read in fhandle:
    if (
      not read.is_unmapped and
      read.query_alignment_length >= 0.8 * ctg_size_map[read.qname]
    ):
      map_ctgs.add(read.qname)
  fhandle.close()
  umap_ctgs = set(ctg_size_map.keys()) - map_ctgs
  return umap_ctgs

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

