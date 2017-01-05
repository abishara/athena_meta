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
    self.logger.log('jointly assemble bins with OLC')

    # collect input contig and local reasm contigs
    self.logger.log('merge input contigs')

    mergedfa_path = os.path.join(self.outdir, 'canu-input-contigs.fa')
    mergedfiltfa_path = os.path.join(self.outdir, 'canu-input-contigs.filt.fa')
    seedsfa_path = os.path.join(self.outdir, 'seed-contigs.fa')

    # load all the bins
    bins = util.load_pickle(self.options.bins_pickle_path)

    #input_paths = [self.options.ctgfasta_path]
    input_paths = []
    seed_ctgs = set()
    for (binid, seed_group) in bins:
      for ctg in seed_group:
        seed_ctgs.add(ctg)
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

    # FIXME uncomment
    util.concat_files(input_paths, mergedfa_path)
    assert is_valid_fasta(mergedfa_path), "merge FASTA not valid"
    #die

    mergedbam_path = os.path.join(self.outdir, 'align-inputs.bam')
    cmd = 'bwa mem -t 8 {} {} | samtools view -bS - | samtools sort -o {} - '.format(
      self.options.ctgfasta_path,
      mergedfa_path,
      mergedbam_path,
    )
    subprocess.check_call(cmd, shell=True)
    cmd = 'samtools index {}'.format(mergedbam_path)
    subprocess.check_call(cmd, shell=True)

    filter_inputs(
      mergedbam_path,
      mergedfa_path,
      mergedfiltfa_path,
    )
    #die


    #self.logger.log('  {} contigs, covering {} bases'.format(
    #  total_asm_contigs,
    #  total_asm_bp,
    #))

    canu0_path = os.path.join(self.outdir, 'canu-asm-1.seeds')
    cmd = \
'{} \
useGrid=1  \
gridOptions="-p owners" \
errorRate=0.06  \
genomeSize=45.00m  \
contigFilter="2 2000 1.0 1.0 2" \
stopOnReadQuality=false \
-assemble \
-d {} \
-p canu \
oeaMemory=12 cnsMemory=32 batMemory=50 \
-pacbio-corrected {}'.format(
      canubin_path,
      canu0_path,
      mergedfiltfa_path
    )
    print 'cmd', cmd
    subprocess.check_call(cmd, shell=True)
    die

    # index assembled contigs 
    self.logger.log('index canu assembled contigs')
    with util.cd(canu0_path):
      #pass
      cmd = 'bwa index canu.contigs.fasta'
      subprocess.check_call(cmd, shell=True)

    # align idba0 contigs to canu contigs
    canu_contigs_path = os.path.join(canu0_path, 'canu.contigs.fasta')

    # align init contigs to canu contigs
    self.logger.log('aligning seed contigs to olc contigs')
    idba0fa_path = self.options.ctgfasta_path
    outsam_path = os.path.join(self.outdir, 'align.on-contig.sam')
    cmd = 'bwa mem -t 4 {} {} > {}'.format(
      canu_contigs_path,
      idba0fa_path,
      outsam_path,
    )
    print 'cmd', cmd
    subprocess.check_call(cmd, shell=True)
    with util.cd(self.outdir):
      print 'cmd', cmd
      cmd = 'cat align.on-contig.sam | samtools view -bS - | samtools sort -o align.on-contig.bam -'
      subprocess.check_call(cmd, shell=True)
      print 'cmd', cmd
      cmd = 'samtools index align.on-contig.bam'
      subprocess.check_call(cmd, shell=True)

    seeds = set() 
    bins = util.load_pickle(self.options.bins_pickle_path)
    for _, _seeds in bins:
      seeds |= set(_seeds)

    idba_fasta = pysam.FastaFile(idba0fa_path)
    self.logger.log('determine unmapped contigs out of {} idba0 contigs'.format(
      idba_fasta.nreferences))

    unmap_idba0_ctgs = get_unmapped_ctgs(
      idba0fa_path,
      os.path.join(self.outdir, 'align.on-contig.bam'))
    self.logger.log('  - {} unmapped'.format(len(unmap_idba0_ctgs)))
    self.logger.log('  - {}/{} seeds unmapped'.format(
      len(unmap_idba0_ctgs & seeds),
      len(seeds),
    ))

    idba0_unmapfa_path = os.path.join(self.outdir, 'idba0.unmap.fa')
    idba0_unmapseedsfa_path = os.path.join(self.outdir, 'idba0.unmap.seeds.fa')


    with open(idba0_unmapfa_path, 'w') as fout1, \
         open(idba0_unmapseedsfa_path, 'w') as fout2:
      for ctg in unmap_idba0_ctgs:
        seq = str(idba_fasta.fetch(ctg).upper())
        fout1.write('>{}\n'.format(ctg))
        fout1.write('{}\n'.format(seq))
        if ctg in seeds:
          fout2.write('>{}\n'.format(ctg))
          fout2.write('{}\n'.format(seq))

    idba0_seedsfa_path = os.path.join(self.outdir, 'idba.seeds.fa')
    with open(idba0_seedsfa_path, 'w') as fout:
      for ctg in seeds:
        seq = str(idba_fasta.fetch(ctg).upper())
        fout.write('>{}\n'.format(ctg))
        fout.write('{}\n'.format(seq))

    final_fa_path = os.path.join(self.outdir, 'final.asm.fa')
    final_fa2_path = os.path.join(self.outdir, 'final.asm2.fa')
    util.concat_files(
      [idba0_unmapfa_path, canu_contigs_path],
      final_fa_path,
    )
    util.concat_files(
      [idba0_unmapseedsfa_path, canu_contigs_path],
      final_fa2_path,
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
    # exclude identity alignments
    if read.qname == fhandle.getrname(read.tid):
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



