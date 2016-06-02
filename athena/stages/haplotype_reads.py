import os
import pysam

#from athena.stages.step import StepChunk
from .step import StepChunk
from ..mlib import util

class HaplotypeReadsStep(StepChunk):
  @staticmethod
  def get_steps(options):
    f = pysam.FastaFile(options.ref_fasta)
    for ctg in f.references:
      l = f.get_reference_length(ctg)
      for i, (b, _) in enumerate(
        util.get_partitions(0, l, options.genome_step_size)
      ):
        if i > 5:
          break
        e = b + options.genome_window_size
        yield HaplotypeReadsStep(options, ctg, b, e)
      break
    f.close()

  def outpaths(self, final=False):
    return {
      "file" : os.path.join(self.output_dir, 'file'),
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
    self.logger.log('i am running')
    file_path = os.path.join(self.output_dir, 'file')
    util.touch(file_path)
    self.logger.log('i am done')

