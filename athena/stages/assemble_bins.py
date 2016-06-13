import os
import pysam
import subprocess
from collections import defaultdict
import glob
import shutil

from ..assembler_tools.architect import architect
from ..assembler_tools.architect import local_assembly

from .step import StepChunk
from ..mlib import util

# FIXME remove harcoded path
idbabin_path = '/home/abishara/sources/idba/bin/idba_ud'

wd = os.path.dirname(os.path.abspath(__file__))
architect_scripts_path = os.path.join(
  wd,
  '..',
  'assembler_tools/architect/scripts/',
)
edgesbin_path = os.path.join(architect_scripts_path, 'pe-connections.py')
containmentbin_path = os.path.join(architect_scripts_path, 'bam_to_containment.py')

class AssembleBinsStep(StepChunk):

  @staticmethod
  def get_steps(options):
    bins = util.load_pickle(options.bins_pickle_path)

    for i, (binid, bcode_set) in enumerate(bins):
      yield AssembleBinsStep(options, binid)

  def outpaths(self, final=False):
    paths = {}
    paths['contig.fa'] = os.path.join(
      self.options.get_bin_dir(self.binid, final=True),
      'contig.fa',
    )
    paths['local-asm-merged.fa'] = os.path.join(
      self.options.get_bin_dir(self.binid, final=True),
      'local-asm-merged.fa',
    )
    paths['architect-scaff.fasta'] = os.path.join(
      self.options.get_bin_dir(self.binid, final=True),
      'scaff.fasta',
    )
    return paths

  @property
  def outdir(self):
    return self.options.get_bin_dir(self.binid, final=True)

  def clean_working(self):
    bin_path = self.options.get_bin_dir(self.binid)
    self.logger.log('removing bin directory {}'.format(bin_path))
    shutil.rmtree(bin_path)
    return 

  def __init__(
    self,
    options,
    binid,
  ):
    self.options = options
    self.binid = binid
    util.mkdir_p(self.outdir)

  def __str__(self):
    ctg, b, e, cidx = self.binid
    return '{}_{}.{}-{}.c{}'.format(
      self.__class__.__name__,
      ctg, b, e, cidx,
    )

  def run(self):
    self.logger.log('assembling barcoded reads for this bin')

    fqdir_path = self.options.get_bin_fq_dir(self.binid)
    allfa_path = os.path.join(fqdir_path, 'tenxreads.fa')
    asmdir_path = self.options.get_bin_asm_dir(self.binid)
    # merge fq fragments to create input reads, create idba *fa reads
    with util.cd(fqdir_path):
      cmd = 'cat *frag.fq > tenxreads.fq'
      subprocess.check_call(cmd, shell=True)
      util.convert_tenx_fq2fa('tenxreads.fq', 'tenxreads.fa')
      util.convert_tenx_fq2bcodefa('tenxreads.fq', 'tenxreads-bcoded.fa')
    # initial idba assembly
    cmd = '{} -r {} -o {}'.format(
      idbabin_path,
      allfa_path,
      asmdir_path,
    )
    subprocess.check_call(cmd, shell=True)
    
    with util.cd(asmdir_path):
      assert os.path.isfile('contig.fa'), 'idba failed to generate contig.fa'

      # index initial contigs
      cmd = 'bwa index contig.fa'
      subprocess.check_call(cmd, shell=True)

      # align reads to contigs
      self.logger.log('aligning reads to contigs')
      cmd = 'bwa mem -p contig.fa {} > align.on-contig.sam'.format(
        '../fqs/tenxreads-bcoded.fa'
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
      cmd = 'samtools index {}'.format(filtbam_path)
      subprocess.check_call(cmd, shell=True)

      # create edges.tsv and containment files from architect
      self.logger.log('generating local reassembly inputs')
      cmd = 'python {} -b {} -f {} -e {}'.format(
        edgesbin_path,
        filtbam_path,
        'contig.fa',
        'edges.tsv',
      )
      subprocess.check_call(cmd, shell=True)
      cmd = 'python {} -b {} -c {}'.format(
        containmentbin_path,
        filtbam_path,
        'containment',
      )
      subprocess.check_call(cmd, shell=True)

      # extract barcoded filtered alignments
      filtfa_path = 'tenxreads-bcoded-filt.fa'
      get_filtered_fa(
        filtbam_path,
        '../fqs/tenxreads-bcoded.fa',
        filtfa_path,
      )

      # run architect for local assemblies
      local_assembly.setup_cwd()
      self.logger.log('performing local reassembly')
      architect.do_local_reassembly(
        'contig.fa',
        'edges.tsv',
        'containment',
        filtfa_path,
      )

      # FIXME for now just run scaffolding in each bucket separately
      architect.do_scaffolding(
        'contig.fa',
        'edges.tsv',
        'containment',
        0,
        'scaff',
      )

    def copyfinal(src, dest):
      shutil.copyfile(
        os.path.join(self.options.get_bin_asm_dir(self.binid), src),
        os.path.join(self.options.get_bin_dir(self.binid, final=True), dest),
      )

    self.logger.log('copying deliverables to final')
    copyfinal('contig.fa', 'contig.fa')
    copyfinal('local-assemblies/local-asm-merged.fa', 'local-asm-merged.fa')
    copyfinal('scaff.fasta', 'scaff.fasta')

    self.logger.log('done')

#--------------------------------------------------------------------------
# helpers
#--------------------------------------------------------------------------
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
    for qname, txt in util.fa_iter(infa_path):
      if qname in qname_set:
        fout.write(txt)

