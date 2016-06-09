import os
import pysam
import subprocess
from ..assembler_tools.haplotyper import haplotyper

#from athena.stages.step import StepChunk
from .step import StepChunk
from ..mlib import util

wd = os.path.dirname(os.path.abspath(__file__))
haplotyperbin_path = os.path.join(
  wd,
  '..',
  'assembler_tools/haplotyper/haplotyper.py',
)
assert os.path.isfile(haplotyperbin_path)

class HaplotypeReadsStep(StepChunk):
  @staticmethod
  def get_steps(options):
    for ctg, b, e in util.get_genome_partitions(
      options.ref_fasta,
      options.genome_window_size,
      options.genome_step_size,
    ):
      if not options.regions or len(options.regions[ctg].find(b,e)) > 0:
        yield HaplotypeReadsStep(options, ctg, b, e)
        break

  def outpaths(self, final=False):
    return {
      'stats' : os.path.join(self.output_dir, 'stats.p'),
    }

  def __init__(
    self,
    options,
    ctg,
    begin,
    end,
  ):
    self.options = options
    self.ctg = ctg
    self.begin = begin
    self.end = end
    self.output_dir = os.path.join(
      self.options.working_dir,
      self.__class__.__name__,
      str(self),
    )
    util.mkdir_p(self.output_dir)

  def __str__(self):
    return '{}_{}.{}-{}'.format(
      self.__class__.__name__,
      self.ctg,
      self.begin,
      self.end,
    )

  def run(self):
    # FIXME make sure to work out of scratch directory
    self.logger.log('clustering reads')
    roi_str = '{}:{}-{}'.format(self.ctg, self.begin, self.end)
    haplotyper.cluster_reads(
      self.options.longranger_bam_path,
      self.options.longranger_vcf_path,
      roi_str,
      self.output_dir,
    )
    self.logger.log('done')
    #file_path = os.path.join(self.output_dir, 'file')
    #util.touch(file_path)

