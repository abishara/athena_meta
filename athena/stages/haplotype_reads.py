import os
import pysam
import subprocess
from ..assembler_tools.haplotyper import haplotyper

from .step import StepChunk
from ..mlib import util

wd = os.path.dirname(os.path.abspath(__file__))
haplotyperbin_path = os.path.join(
  wd,
  '..',
  'assembler_tools/haplotyper/haplotyper.py',
)
assert os.path.isfile(haplotyperbin_path)

MIN_CONTIG_SIZE = 20000

#--------------------------------------------------------------------------
# base class
#--------------------------------------------------------------------------
class BaseStep(StepChunk):

  @classmethod
  def get_steps(cls, options):

    ctg_size_map = util.get_fasta_sizes(options.ctgfasta_path)
    if os.path.isfile(options.bins_pickle_path):
      bins = util.load_pickle(options.bins_pickle_path)
    else:
      print 'creating haplotyping bins'
      bins = map(
        lambda(b, v): ('bin.{}'.format(b), v),
        enumerate(util.get_fasta_partitions(options.ctgfasta_path, 128)),
      )
      print '  - saving {}'.format(options.bins_pickle_path)
      util.write_pickle(options.bins_pickle_path, bins)

    total_ctgs = 0
    total_bases = 0
    for binid, group in bins:
      total_bases += sum(map(lambda(c): ctg_size_map[c], group))
      total_ctgs += len(group)
      yield cls(options, binid, group)

    print 'grouped {} total_ctgs covering {} bases'.format(
      total_ctgs,
      total_bases)

  @property
  def outdir(self):
    return os.path.join(
      self.options.results_dir,
      self.__class__.__name__,
      str(self),
    )
 
  def __str__(self):
    return '{}_{}'.format(
      self.__class__.__name__,
      self.binid,
    )

  def __init__(
    self,
    options,
    binid,
    ctgs,
  ):
    self.options = options
    self.ctgs = ctgs
    self.binid = binid
    util.mkdir_p(self.outdir)
    util.mkdir_p(self.options.get_bin_dir(self.binid))

#--------------------------------------------------------------------------
# call variants
#--------------------------------------------------------------------------
# FIXME hardcoded path
freebayesbin_path = '/home/abishara/sources/freebayes/bin/freebayes'

class CallVariantsStep(BaseStep):

  def outpaths(self, final=False):
    paths = {}
    paths['vcf'] = os.path.join(self.outdir, 'vars.vcf.gz')
    paths['vcf_index'] = os.path.join(self.outdir, 'vars.vcf.gz.tbi')
    return paths

  def run(self):
    self.logger.log('creating regions bed file')
    ctg_size_map = util.get_fasta_sizes(self.options.ctgfasta_path)
    bindir_path = self.options.get_bin_dir(self.binid)
    bed_path = os.path.join(bindir_path, 'roi.bed')
    with open(bed_path, 'w') as fout:
      for ctg in self.ctgs:
        fout.write('{}\t{}\t{}\n'.format(ctg, 0, ctg_size_map[ctg]))
      
    self.logger.log('freebayes for contigs in this group')
    vcf_path = os.path.join(bindir_path, 'vars.filt.vcf')
    with open(vcf_path, 'w') as fout:
      #-Y 100 -C 4 -F 0.1 
      subprocess.check_call([
        freebayesbin_path,
        '-f', self.options.ctgfasta_path,
        '-b', self.options.reads_ctg_bam_path,
        '-0',
        #'-Y', '120', '-C', '4', '-F', '0.1', '-q', '20',
        '-t', bed_path], stdout=fout)

    self.logger.log('bgzip and tabix index')
    vcfgz_path = os.path.join(self.outdir, 'vars.vcf.gz')
    cmd = 'bgzip -c {} > {}'.format(vcf_path, vcfgz_path)
    subprocess.check_call(cmd, shell=True)
    cmd = 'tabix -p vcf {}'.format(vcfgz_path)
    subprocess.check_call(cmd, shell=True)

    self.logger.log('done')


#--------------------------------------------------------------------------
# haplotype reads
#--------------------------------------------------------------------------
class HaplotypeReadsStep(BaseStep):

  def outpaths(self, final=False):
    paths = {}
    paths['stats'] = os.path.join(self.outdir, 'stats.p')
    paths['strains'] = os.path.join(self.outdir, 'strains.vcf')
    return paths

  def run(self):
    var_step = CallVariantsStep(self.options, self.binid, self.ctgs)
    invcf_path = var_step.outpaths()['vcf']
    self.logger.log('haplotyping contigs in group')

    stats_paths = []
    hapvcf_paths = []
    ctg_size_map = util.get_fasta_sizes(self.options.ctgfasta_path)
    for ctg in self.ctgs:
      self.logger.log('  - {}'.format(ctg))

      ## FIXME remove
      #roi_str = '{}:{}-{}'.format(ctg, 0, 1000)
      roi_str = '{}:{}-{}'.format(ctg, 0, ctg_size_map[ctg])
      outdir = os.path.join(self.outdir, ctg)
      haplotyper.cluster_reads(
        self.options.reads_ctg_bam_path,
        invcf_path,
        roi_str,
        outdir,
      )
      stats_path = os.path.join(outdir, 'stats.p')
      hapvcf_path = os.path.join(outdir, 'clusters.vcf')
      stats_paths.append((ctg, stats_path))
      hapvcf_paths.append((ctg, hapvcf_path))

    self.logger.log('merging outputs from all contigs in groups')
    strains_vcf_path = os.path.join(self.outdir, 'strains.vcf')
    with open(strains_vcf_path, 'w') as fout:
      header = False
      for _, vcf_path in hapvcf_paths:
        with open(vcf_path) as fin:
          for line in fin:
            if line.startswith('#'):
              if not header:
                fout.write(line)
              continue
            fout.write(line)
        header = True
    
    mstats = [] 
    merged_stats_path = os.path.join(self.outdir, 'stats.p')
    for ctg, stats_path in stats_paths:
      stats = util.load_pickle(stats_path)
      mstats.append((ctg, stats))
    util.write_pickle(merged_stats_path, mstats)

    self.logger.log('done')

#--------------------------------------------------------------------------
# helpers
#--------------------------------------------------------------------------
def get_vcf_path(vcf_path, ctg):
  vcf_fhandle = vcf.Reader(filename=vcf_path)
  try:
    vcf_iter = vcf_fhandle.fetch(ctg)
  except Exception as e:
    print 'exception parsing vcf', str(e)
    return []
  return vcf_iter

