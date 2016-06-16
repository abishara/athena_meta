
import os
import pysam
import shutil

from .step import StepChunk
from ..mlib import util
import assemble_bins

#--------------------------------------------------------------------------
# base class
#--------------------------------------------------------------------------
class AssembleGroupsStep(StepChunk):

  @staticmethod
  def get_steps(options):
    groups = util.load_pickle(options.groups2_pickle_path)

    for gid, bins in groups:
      yield AssembleGroupsStep(options, gid, bins)

  @property
  def outdir(self):
    return self.options.get_group_dir(self.gid, final=True)

  def outpaths(self, final=False):
    paths = {}
    paths['olc-asm.contig.fasta'] = os.path.join(self.outdir, 'olc-asm.contig.fasta')
    paths['olc-asm.scaff.fasta'] = os.path.join(self.outdir, 'olc-asm.scaff.fasta')
    return paths

  def clean_working(self):
    group_path = self.options.get_group_dir(gid)
    self.logger.log('removing group directory {}'.format(group_path))
    shutil.rmtree(group_path)
    return 

  def __init__(
    self,
    options,
    gid,
    bins,
  ):
    self.options = options
    self.gid = gid
    self.bins = bins
    util.mkdir_p(self.outdir)

  def __str__(self):
    return '{}.group.{}'.format(self.__class__.__name__, self.gid)

  def run(self):
    self.logger.log('jointly assemble bins with OLC and architect scaffolding')

    # collect input contig and local reasm contigs
    for (binid, bcode_set) in self.bins:
      asm_step = assemble_bins.AssembleBinnedStep(self.options, binid)
      contig_path = asm_step.outpaths()['contig.fa']
      localasm_path = asm_step.outpaths()['local-asm-merged.fa']
      self.logger.log('for binid {}'.format(binid))
      self.logger.log('  - reading contig {}'.format(contig_path))
      self.logger.log('  - reading local-asm {}'.format(localasm_path))

    self.logger.log('lol assembling')
    self.logger.log('done')

