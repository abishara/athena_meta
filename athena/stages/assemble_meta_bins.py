import os
import pysam
import subprocess
from collections import defaultdict
import glob
import shutil
from itertools import izip
import random

from ..assembler_tools.architect import architect
from ..assembler_tools.architect import local_assembly

from .step import StepChunk
from ..mlib import util

# NOTE must be in path
idbabin_path = 'idba_ud'

wd = os.path.dirname(os.path.abspath(__file__))
architect_scripts_path = os.path.join(
  wd,
  '..',
  'assembler_tools/architect/scripts/',
)
edgesbin_path = os.path.join(architect_scripts_path, 'pe-connections.py')
containmentbin_path = os.path.join(architect_scripts_path, 'bam_to_containment.py')

class AssembleMetaBinnedStep(StepChunk):

  @staticmethod
  def get_steps(options):
    bins = util.load_pickle(options.bins_pickle_path)
    bins2 = util.load_pickle(options.bins2_pickle_path)
    assert len(bins) == len(bins2)
    
    for i, ((binid, bcode_set), (_, ds_bcode_set, hint_ctgs)) in \
        enumerate(izip(bins, bins2)):
      yield AssembleMetaBinnedStep(options, binid, bcode_set, hint_ctgs, ds_bcode_set)

  def __init__(
    self,
    options,
    binid,
    bcode_set,
    hint_ctgs,
    ds_bcode_set=None,
  ):
    self.options = options
    self.binid = binid
    self.bcode_set = bcode_set
    self.hint_ctgs = hint_ctgs
    self.ds_bcode_set = ds_bcode_set
    util.mkdir_p(self.outdir)

  def __str__(self):
    return '{}.{}'.format(self.__class__.__name__, self.binid)

  @property
  def outdir(self):
    return self.options.get_bin_dir(self.binid, final=True)

  def outpaths(self, final=False):
    paths = {}
    #paths['contig.fa'] = os.path.join(self.outdir, 'contig.fa')
    #paths['shit'] = 'shit'
    paths['local-asm-merged.fa'] = os.path.join(self.outdir, 'local-asm-merged.fa')
    return paths

  def clean_working(self):
    bin_path = self.options.get_bin_dir(self.binid)
    self.logger.log('removing bin directory {}'.format(bin_path))
    shutil.rmtree(bin_path)
    return 

  def run(self):

    self.logger.log('assembling barcoded reads for this bin')
    self.logger.log('  - {} seed contigs'.format(len(self.hint_ctgs)))
    print self.hint_ctgs

    fqdir_path = self.options.get_bin_fq_dir(self.binid)
    # merge fragments
    self.logger.log('merging input fq fragments')
    with util.cd(fqdir_path):
      cmd = 'cat *frag.fq > tenxreads.fq'
      subprocess.check_call(cmd, shell=True)

    allfa_path = os.path.join(fqdir_path, 'tenxreads-bcoded.fa')
    downsamplefa_path = os.path.join(fqdir_path, 'tenxreads-bcoded.ds.fa')
    lrhintsfa_path = os.path.join(fqdir_path, 'lr-hints.fa')
    asmdir_path = self.options.get_bin_asm_dir(self.binid)
    with util.cd(fqdir_path):
      # create barcoded idba *fa reads
      util.convert_tenx_fq2bcodefa('tenxreads.fq', 'tenxreads-bcoded.fa')

      ## downsample for faster idba assembly if needed
      #if self.ds_bcode_set != None:
      #  self.logger.log('downsample reads')
      #  self.logger.log('  - orig barcodes {}, new barcodes {}'.format(
      #    len(self.bcode_set),
      #    len(self.ds_bcode_set),
      #  ))
      #  with open('tenxreads-bcoded.ds.fa', 'w') as fout:
      #    for e in util.fa_iter('tenxreads-bcoded.fa'):
      #      bcode = e.qname.split('$')[0]
      #      if bcode in self.ds_bcode_set:
      #        fout.write(e.txt)
      #else:
      #  pass
      #  shutil.copyfile('tenxreads-bcoded.fa', 'tenxreads-bcoded.ds.fa')

      # create long read hints
      with open('lr-hints.fa', 'w') as fout, \
           open('lr-hints2.fa', 'w') as fout2:
        ctg_fasta = pysam.FastaFile(self.options.ctgfasta_path)
        for ctg in self.hint_ctgs:
          seq = str(ctg_fasta.fetch(ctg).upper())
          fout2.write('>{}\n{}\n'.format(ctg, seq))
          for _ in xrange(10):
            fout.write('>{}\n{}\n'.format(ctg, seq))

    self.logger.log('performing idba assembly')

    # initial idba assembly
    #self.do_idba_assemble(
    #  downsamplefa_path,
    #  asmdir_path,
    #  lrhintsfa_path,
    #)
    util.mkdir_p(asmdir_path)

    with util.cd(asmdir_path):

      shutil.copy('../fqs/lr-hints2.fa', '.')
      cmd = 'bwa index lr-hints2.fa'
      subprocess.check_call(cmd, shell=True)

      # align reads to contigs
      self.logger.log('aligning reads to contigs')
      cmd = 'bwa mem -p lr-hints2.fa {} > align.on-contig.sam'.format(
        '../fqs/tenxreads-bcoded.fa'
      )
      subprocess.check_call(cmd, shell=True)
      cmd = 'cat align.on-contig.sam | samtools view -bS - | samtools sort - align.on-contig.sorted'
      subprocess.check_call(cmd, shell=True)

      # filter bad reads, create idba *fa reads with only filtered reads
      allbam_path = 'align.on-contig.sorted.bam'
      #filtbam_path = 'align.on-contig.sorted.filt.bam'
      #filter_alignments(
      #  allbam_path,
      #  filtbam_path,
      #)
      cmd = 'samtools index {}'.format(allbam_path)
      subprocess.check_call(cmd, shell=True)

      # create edges.tsv and containment files from architect
      self.logger.log('generating local reassembly inputs')
      cmd = 'python {} -b {} -f {} -e {}'.format(
        edgesbin_path,
        allbam_path,
        'lr-hints2.fa',
        'edges.tsv',
      )
      subprocess.check_call(cmd, shell=True)
      cmd = 'python {} -b {} -c {}'.format(
        containmentbin_path,
        allbam_path,
        'containment',
      )
      subprocess.check_call(cmd, shell=True)

      # run architect for local assemblies
      local_assembly.setup_cwd()
      self.logger.log('performing local reassembly')
      architect.do_local_reassembly(
        'lr-hints2.fa',
        'edges.tsv',
        'containment',
        '../fqs/tenxreads-bcoded.fa',
      )
    
    def copyfinal(src, dest):
      shutil.copyfile(
        os.path.join(self.options.get_bin_asm_dir(self.binid), src),
        os.path.join(self.options.get_bin_dir(self.binid, final=True), dest),
      )

    self.logger.log('copying deliverables to final')
    copyfinal('local-assemblies/local-asm-merged.fa', 'local-asm-merged.fa')

    self.logger.log('done')

#--------------------------------------------------------------------------
# helpers
#--------------------------------------------------------------------------
  def do_idba_assemble(
    self,
    downsamplefa_path,
    asmdir_path,
    lrhintsfa_path=None,
  ):
    cmd = '{} -r {} -l {} -o {}'.format(
      idbabin_path,
      downsamplefa_path,
      lrhintsfa_path,
      asmdir_path,
    )
    cmd = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
    retcode = cmd.wait()
    print "::::", retcode
    contigs_path = os.path.join(asmdir_path, 'contig.fa')
    if not os.path.exists(contigs_path):
      if retcode != 0 and  "invalid insert distance" in cmd.stderr.read():
        self.logger.log("idba_ud internal error; this is probably a bug in idba_ud, so we have to bail here")
        open(contigs_path, "w")
      else:
        error_message = "assembly failed to produce contig.fa"
        self.logger.log(error_message)
        raise Exception()
    elif retcode != 0:
      self.logger.log(
        "something funky happened while running idba_ud (got error code "
        "{}) but there's not much we can do so continuing".format(retcode))

def filter_alignments(in_path, out_path):
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

def get_filtered_fa(
  filtbam_path,
  infa_path,
  outfa_path,
):
  # load all filtered query names to accept
  qname_set = set()
  fin = pysam.Samfile(filtbam_path, 'rb')
  for read in fin:
    qname_set.add(read.qname)
  fin.close()

  with open(outfa_path, 'w') as fout:
    for e in util.fa_iter(infa_path):
      if e.qname_full in qname_set:
        fout.write(e.txt)

  
