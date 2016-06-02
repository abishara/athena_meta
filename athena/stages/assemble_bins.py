import os
import pysam
import subprocess
from collections import defaultdict
import glob

#from athena.stages.step import StepChunk
from .step import StepChunk
from ..mlib import util

# FIXME remove hardcoded path
idbabin_path = '/home/abishara/projs/idba-rc/bin/idba_ud'

class AssembleBinsStep(StepChunk):

  @staticmethod
  def get_steps(options):
    bins = util.load_pickle(self.options.bins_pickle_path)

    for i, (binid, bcode_set) in enumerate(bins):
      yield AssembleBinsStep(options, binid)

  def outpaths(self, final=False):
    return {
    }

  def __init__(
    self,
    options,
    binid,
  ):
    self.options = options
    self.binid = binid

  def __str__(self):
    ctg, b, e, cidx = binid
    return '{}_{}.{}-{}.c{}'.format(
      self.__class__.__name__,
      ctg, b, e, cidx,
    )

  def run(self):
    # FIXME make sure to work out of scratch directory
    self.logger.log('assembling barcoded reads for this bin')

    fqdir_path = self.options.get_bin_fq_dir(binid)
    allfa_path = os.path.join(fqdir_path, 'tenxreads.fa')
    asmdir_path = self.options.get_bin_asm_dir(binid)
    # merge fq fragments to create input reads, create idba *fa reads
    with util.cd(fqdir_path):
      cmd = 'cat *frag.fq > tenxreads.fq'
      subprocess.check_call(cmd, shell=True)
      util.convert_tenx_fq2fa('tenxreads.fq', 'tenxreads.fa')

    # initial idba assembly
    cmd = '{} -r {} -o {}'.format(
      idbabin_path,
      allfa_path,
      asmdir_path,
    )
    subprocess.check_call(cmd, shell=True)
    
    with util.cd(asmdir_path):
      assert os.path.isfile(contig.fa), 'idba failed to generate contig.fa'

      # index initial contigs
      cmd = 'bwa index contig.fa'
      subprocess.check_call(cmd, shell=True)

      # align reads to contigs
      cmd = 'bwa mem -p contig.fa {} > align.on-contig.sam'.format(
        '../fqs/tenxreads.fa'
      )
      subprocess.check_call(cmd, shell=True)
      cmd = 'cat align.on-contig.sam | samtools view -bS - | samtools sort - align.on-contig.sorted'
      subprocess.check_call(cmd, shell=True)

      # filter bad reads, create idba *fa reads with only filtered reads
      allbam_path = 'align.on-contig.sorted.bam'
      filtbam_path = 'align.on-contig.sorted.filt.bam'
      filter_alignments(
        allbam_path,
        filtbam_path,
      )

      # create edges.tsv and containment files from architect
      # run architect for local assemblies

    self.logger.log('done')

#--------------------------------------------------------------------------
# helpers
#--------------------------------------------------------------------------
def filter_alignments(inbam_path, outbam_path):
  def get_clip_info(read):
    (c, l) = read.cigar[0]
    lclip = (c in [4,5])
    (c, l) = read.cigar[-1]
    rclip = (c in [4,5])
  
    return (lclip, rclip)
  
  fin = pysam.Samfile(in_path, 'rb')
  fullrid_set = set()
  for read in fin:
    rid = (read.qname, read.is_read1)
    if read.is_unmapped:
      continue
    (lclip, rclip) = get_clip_info(read)
    if not (lclip or rclip):
      fullrid_set.add(rid)
  fin.close()
  
  fin = pysam.Samfile(in_path, 'rb')
  fout = pysam.Samfile(out_path, 'wb', template=fin)
  for read in fin:
    rid = (read.qname, read.is_read1)
    prid = (read.qname, not read.is_read1)
    if rid in fullrid_set or prid in fullrid_set:
      fout.write(read)
    #if rid in fullrid_set and prid in fullrid_set:
    #  fout.write(read)
  fin.close()
  fout.close()
  return


