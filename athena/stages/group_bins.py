import os
import pysam
import shutil

from .step import StepChunk
from ..mlib import util

class GroupBinsStep(StepChunk):

  @staticmethod
  def get_steps(options):
    yield GroupBinsStep(options)

  def __init__(
    self,
    options,
  ):
    self.options = options
    util.mkdir_p(self.outdir)

  def __str__(self):
    return self.__class__.__name__

  @property
  def outdir(self):
    return os.path.dirname(self.options.groups_pickle_path)

  def outpaths(self, final=False):
    paths = {}
    paths['groups.p'] = self.options.groups_pickle_path
    return paths

  def run(self):
    self.logger.log('determine bins to jointly assemble')

    bins = util.load_pickle(self.options.bins_pickle_path)

    # FIXME for now just put each bin into its own group
    groups = []
    groups2 = []
    for gid, (binid, bcode_set) in enumerate(bins):
      groups.append((gid, bcode_set))
      groups2.append((gid, [(binid, bcode_set)]))

    util.write_pickle(self.options.groups_pickle_path, groups)
    util.write_pickle(self.options.groups2_pickle_path, groups2)

    self.logger.log('done')

