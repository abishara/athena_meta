import os
import pysam
import subprocess
from collections import defaultdict
import glob

from ..assembler_tools.haplotyper import haplotyper

from .step import StepChunk
from ..mlib import util
import haplotype_reads

class CollectReadsStep(StepChunk):

  @staticmethod
  def get_steps(options):
    # determine reads to collect from haplotyper step if not generated
    if not os.path.isfile(options.bins_pickle_path):
      bins = []

      for hap_step in haplotype_reads.HaplotypeReadsStep.get_steps(options):
        hapout_path = hap_step.outpaths()['stats']
        (cluster_info_map, _) = util.load_pickle(hapout_path)
        for cidx, info in cluster_info_map.items():
          numreads, bcode_set, _, assm = info
          ctg, b, e = hap_step.ctg, hap_step.begin, hap_step.end
          binid = (ctg, b, e, cidx)
          if assm and cidx != None:
            bins.append((binid, bcode_set))
            # ensure output fq frag directory exists
            util.mkdir_p(options.get_bin_fq_dir(binid))
      util.write_pickle(options.bins_pickle_path, bins)

    # strip over fastqs to load all fq fragments
    rootfq_path = options.longranger_fqs_path
    for fq_path in glob.glob(rootfq_path + '/chnk*/files/*fastq*gz'):
      yield CollectReadsStep(options, fq_path)

  def outpaths(self, final=False):
    paths = {}
    paths['pass.file'] = os.path.join(self.outdir, 'pass')
    return paths
 
  @property
  def outdir(self):
    return os.path.join(
      self.options.results_dir,
      self.__class__.__name__,
      str(self),
    )

  def __init__(
    self,
    options,
    fq_path,
  ):
    self.options = options
    self.fq_path = fq_path
    util.mkdir_p(self.outdir)

  def __fqid(self):
    return os.path.basename(os.path.dirname(os.path.dirname(self.fq_path)))

  def __str__(self):
    return '{}_{}'.format(
      self.__class__.__name__,
      self.__fqid(),
    )

  def run(self):
    self.logger.log('collecting reads for each bin')

    allbcode_set = set()
    bcode_groups_map = defaultdict(set)
    groupf_map = {}
    bins = util.load_pickle(self.options.bins_pickle_path)

    # open a file handle for each bin
    for binid, bcode_set in bins:
      fqfrag_path = os.path.join(
        self.options.get_bin_fq_dir(binid),
        '{}.frag.fq'.format(self.__fqid()),
      )
      groupf_map[binid] = open(fqfrag_path, 'w')
      allbcode_set |= bcode_set
      for bcode in bcode_set:
        bcode_groups_map[bcode].add(binid)

    # for each barcoded read, write to all bins that have that barcode
    for bcode, rtxt in util.tenx_fastq_iter(self.fq_path):
      if bcode in allbcode_set:
        for binid in bcode_groups_map[bcode]:
          groupf_map[binid].write(rtxt)

    for f in groupf_map.values():
      f.close()

    passfile_path = os.path.join(self.outdir, 'pass')
    util.touch(passfile_path)
    self.logger.log('done')

