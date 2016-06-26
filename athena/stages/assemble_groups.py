import os
import pysam
import subprocess
import shutil
from collections import defaultdict, Counter

from .step import StepChunk
from ..mlib import util
import assemble_bins

# FIXME harcoded paths
# NOTE must be in path

canubin_path = '/home/abishara/sources/canu/Linux-amd64/bin/canu'
olapbin_path = '/home/abishara/sources/canu/Linux-amd64/bin/ovStoreDump'

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

    asmdir_path = self.options.get_group_asm_dir(self.gid)
    util.mkdir_p(asmdir_path)

    # collect input contig and local reasm contigs
    self.logger.log('merge input contigs')
    mergedfa_path = os.path.join(asmdir_path, 'input-contigs.fa')
    origmergedfa_path = os.path.join(asmdir_path, 'orig-input-contigs.fa')
    mergedseeds_path = os.path.join(asmdir_path, 'seeds.fa')
    mergedlocal_path = os.path.join(asmdir_path, 'local.fa')
    paths = []
    contig_paths = []
    local_paths = []
    for (binid, bcode_set) in self.bins:
      asm_step = assemble_bins.AssembleBinnedStep(self.options, binid)
      paths.append(asm_step.outpaths()['contig.fa'])
      paths.append(asm_step.outpaths()['local-asm-merged.fa'])
      contig_paths.append(asm_step.outpaths()['contig.fa'])
      local_paths.append(asm_step.outpaths()['local-asm-merged.fa'])
    util.concat_files(paths, mergedfa_path)
    util.concat_files(paths, origmergedfa_path)
    util.concat_files(contig_paths, mergedseeds_path)
    util.concat_files(local_paths, mergedlocal_path)

    with util.cd(asmdir_path):
      canu0_path = 'canu-asm-0'
      canu1_path = 'canu-asm-1'
      mergedfa_path = 'input-contigs.fa'

      # align reads to contigs
      # FIXME remove
      # index initial contigs
      cmd = 'bwa index seeds.fa'
      #subprocess.check_call(cmd, shell=True)
      self.logger.log('aligning reads to contigs')
      cmd = 'bwa mem seeds.fa local.fa > align.on-seeds.sam'
      #subprocess.check_call(cmd, shell=True)
      cmd = 'cat align.on-seeds.sam | samtools view -bS - | samtools sort - align.on-seeds.sorted'
      #subprocess.check_call(cmd, shell=True)
      filtqname_set = get_hq_locals('align.on-seeds.sorted.bam')

      get_filtered_fa(
        filtqname_set,
        'local.fa',
        'local.filt.fa',
      )
      # FIXME remove end
      util.concat_files(['seeds.fa', 'local.filt.fa'], mergedfa_path)

      self.logger.log('initial canu run')
      cmd = '{} useGrid=0 errorRate=0.01 genomeSize=0.02m stopOnReadQuality=false -d {} -p canu -pacbio-corrected {}'.format(
        canubin_path,
        canu0_path,
        mergedfa_path,
      )
      #print 'cmd', cmd
      #subprocess.check_call(cmd, shell=True)

      self.logger.log('assemble reads trimmed by {500, 1000}bp')
      trim500fa_path = 'input-contigs.trim500.fa'
      trim1000fa_path = 'input-contigs.trim1000.fa'
      origtrim1000fa_path = 'input-contigs.orig.trim1000.fa'
      origtrim2000fa_path = 'input-contigs.orig.trim2000.fa'
      trim_reads({},{},
        mergedfa_path,
        trim500fa_path, trim_all=True, trim_size=500)
      trim_reads({},{},
        mergedfa_path,
        trim1000fa_path, trim_all=True, trim_size=1000)
      trim_reads({},{},
        origmergedfa_path,
        origtrim1000fa_path, trim_all=True, trim_size=1000)
      trim_reads({},{},
        origmergedfa_path,
        origtrim2000fa_path, trim_all=True, trim_size=2000)

      cmd = '{} useGrid=0 errorRate=0.01 genomeSize=0.02m stopOnReadQuality=false -d {} -p canu -pacbio-corrected {}'.format(
        canubin_path,
        'canu-alex-asm',
        'alex-asm/alex-asm.inputs.fa',
      )
      print 'cmd', cmd
      subprocess.check_call(cmd, shell=True)
      die
      cmd = '{} useGrid=0 errorRate=0.01 genomeSize=0.02m stopOnReadQuality=false -d {} -p canu -pacbio-corrected {}'.format(
        canubin_path,
        'canu-asm-trim2000-orig',
        origtrim2000fa_path,
      )
      print 'cmd', cmd
      subprocess.check_call(cmd, shell=True)
      cmd = '{} useGrid=0 errorRate=0.01 genomeSize=0.02m stopOnReadQuality=false -d {} -p canu -pacbio-corrected {}'.format(
        canubin_path,
        'canu-asm-trim1000-orig',
        origtrim1000fa_path,
      )
      print 'cmd', cmd
      subprocess.check_call(cmd, shell=True)
      cmd = '{} useGrid=0 errorRate=0.01 genomeSize=0.02m stopOnReadQuality=false -d {} -p canu -pacbio-corrected {}'.format(
        canubin_path,
        'canu-asm-trim500',
        trim500fa_path,
      )
      print 'cmd', cmd
      subprocess.check_call(cmd, shell=True)
      cmd = '{} useGrid=0 errorRate=0.01 genomeSize=0.02m stopOnReadQuality=false -d {} -p canu -pacbio-corrected {}'.format(
        canubin_path,
        'canu-asm-trim1000',
        trim1000fa_path,
      )
      print 'cmd', cmd
      subprocess.check_call(cmd, shell=True)

      die









      self.logger.log('selectively trim reads')
      trimmedfa_path = 'input-contigs.trimmed.fa'
      (qname_3p_overlaps, qname_5p_overlaps) = get_overlaps(canu0_path)
      pre_dolap_qname_set = set(qname_3p_overlaps.keys()) & set(qname_5p_overlaps.keys())
      pre_solap_qname_set = set(qname_3p_overlaps.keys()) | set(qname_5p_overlaps.keys())
      self.logger.log('num dolaps {} solaps {}'.format(
        len(pre_dolap_qname_set),
        len(pre_solap_qname_set),
      ))
      trim_reads(
        qname_3p_overlaps,
        qname_5p_overlaps,
        mergedfa_path,
        trimmedfa_path,
      )

      self.logger.log('secondary canu run')
      cmd = '{} useGrid=0 errorRate=0.01 genomeSize=0.02m stopOnReadQuality=false -d {} -p canu -pacbio-corrected {}'.format(
        canubin_path,
        canu1_path,
        trimmedfa_path,
      )
      print 'cmd', cmd
      subprocess.check_call(cmd, shell=True)

      # dump trimmed reads that have dual and single overlaps
      (qname_3p_overlaps, qname_5p_overlaps) = get_overlaps(canu1_path)
      dolap_qname_set = set(qname_3p_overlaps.keys()) & set(qname_5p_overlaps.keys())
      solap_qname_set = set(qname_3p_overlaps.keys()) | set(qname_5p_overlaps.keys())
      self.logger.log('post num dolaps {} solaps {}'.format(
        len(dolap_qname_set),
        len(solap_qname_set),
      ))
      print 'removed from dolaps', (pre_dolap_qname_set - dolap_qname_set)
      print 'removed from solaps', (pre_solap_qname_set - solap_qname_set)


    self.logger.log('done')

#--------------------------------------------------------------------------
# helpers
#--------------------------------------------------------------------------
def trim_reads(
  qname_3p_overlaps,
  qname_5p_overlaps,
  infa_path,
  outfa_path,
  trim_all=False,
  trim_size=500,
):
  def trim3p(seq): return seq[:-trim_size]
  def trim5p(seq): return seq[trim_size:]

  with open(outfa_path, 'w') as fout:
    for e in util.fa_iter(infa_path):
      seq = e.seq
      header_aux = ''
      if trim_all or (e.qname not in qname_3p_overlaps):
        seq = trim3p(seq)
        header_aux += ' trim3p'
      if trim_all or (e.qname not in qname_5p_overlaps):
        seq = trim5p(seq)
        header_aux += ' trim5p'
      fout.write(e.header+header_aux+'\n')
      fout.write(seq+'\n')

def get_filtered_fa(
  qname_set,
  infa_path,
  outfa_path,
):
  with open(outfa_path, 'w') as fout:
    for e in util.fa_iter(infa_path):
      if e.qname in qname_set:
        fout.write(e.txt)

def get_overlaps(canu_path):
  readnames_path = os.path.join(canu_path, 'unitigging/canu.gkpStore/readNames.txt')
  edges_path = os.path.join(canu_path, 'unitigging/4-unitigger/canu.best.edges')
  rid_map = {}
  with open(readnames_path) as f:
    for line in f:
      rid, rest = line.split('\t')
      qname = rest.split()[0]
      rid_map[rid] = qname

  # invoke overlaps for all reads
  gkp_path = os.path.join(canu_path, 'unitigging/canu.gkpStore')
  ovl_path = os.path.join(canu_path, 'unitigging/canu.ovlStore')
  qname_3p_overlaps = defaultdict(None)
  qname_5p_overlaps = defaultdict(None)
  util.mkdir_p('olaps')
  for rid in rid_map:
    outty = os.path.join('olaps', '{}.olap.txt'.format(rid_map[rid]))
    outty2 = os.path.join('olaps', '{}.olap2.txt'.format(rid_map[rid]))
    with open(outty, 'w') as stdoutf:
      cmd = '{} -G {} -O {} -p {}'.format(
        olapbin_path,
        gkp_path,
        ovl_path,
        rid,
      )
      pp = subprocess.Popen(cmd, shell=True, stdout=stdoutf)
      pp.wait()

    with open(outty) as fin, \
         open(outty2, 'w') as fout:
      first = True
      rb, re = None, None
      for line in fin.readlines():
        words = line.split()
        hrid = words[0]
        hqname = rid_map[hrid]
        newwords = ['{0:40s}'.format(hqname)]
        newwords.extend(words[1:])
        newline = ' '.join(newwords)
        fout.write(newline+'\n')
        if first:
          first = False
          assert rid == hrid
          rb = int(words[2])
          re = int(words[3])
        else:
          mb = int(words[2])
          me = int(words[3])
          if mb == rb:
            qname_5p_overlaps[rid_map[rid]] = rid_map[hrid]
          if me == re:
            qname_3p_overlaps[rid_map[rid]] = rid_map[hrid]

  return (
    qname_3p_overlaps,
    qname_5p_overlaps,
  )

  #qname_3p_overlaps = defaultdict(None)
  #qname_5p_overlaps = defaultdict(None)
  #with open(edges_path) as f:
  #  for line in f:
  #    if line.startswith('#'):
  #      continue
  #    words = line.split('\t')
  #    rid = words[0]
  #    p5hit = words[2]
  #    p3hit = words[4]
  #    contained = (words[-1] == 'contained')
  #    qname = rid_map[rid]
  #    if p5hit != '0':
  #      qname_5p_overlaps[qname] = rid_map[p5hit]
  #    if p3hit != '0':
  #      qname_3p_overlaps[qname] = rid_map[p3hit]
  #return (
  #  qname_3p_overlaps,
  #  qname_5p_overlaps,
  #)


def get_hq_locals(in_path):
  qname_count = Counter()
  qname_overhang_set = set()
  fin = pysam.Samfile(in_path, 'rb')
  for read in fin:
    if read.is_unmapped:
      continue
    if read.query_alignment_length < 500:
      continue
    qname_count[read.qname] += 1
    eslop = read.query_length - read.query_alignment_end
    bslop = read.query_alignment_start
    if eslop >= 750 or bslop >= 750:
      qname_overhang_set.add(read.qname)

  hq_qnames = set()
  for qname, cnt in qname_count.items():
    if cnt > 1 or qname in qname_overhang_set:
      hq_qnames.add(qname)

  fin.close()
  return hq_qnames

