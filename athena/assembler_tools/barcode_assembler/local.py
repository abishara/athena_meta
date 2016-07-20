import os
import pysam
import subprocess
from collections import defaultdict, Counter
import glob
import random
import numpy as np
from bx.intervals.cluster import ClusterTree

from athena.mlib import util
from athena.mlib.fq_idx import FastqIndex

# NOTE must be in path
idbabin_path = 'idba_ud'

MIN_SEED_SIZE = 500
MIN_SEED_SIZE = 400

SEED_SELF_ASM_SIZE = 10000

class LocalAssembly(object):

  def __init__(
    self,
    uid,
    root_ctg,
    root_pos,
    link_ctg,
    link_pos,
    bcode_set,
    ds_bcode_set,
  ):
    self.uid = uid
    self.root_ctg = root_ctg
    self.root_pos = root_pos
    self.link_ctg = link_ctg
    self.link_pos = link_pos
    self.bcode_set = bcode_set
    self.ds_bcode_set = ds_bcode_set

  def __str__(self):
    return '{}_{}.{}'.format(self.root_ctg, self.link_ctg, self.uid)

class LocalAssembler(object):

  def __init__(
    self,
    root_ctg,
    ctgfasta_path,
    reads_ctg_bam_path,
    longranger_fqs_path,
    asmrootdir_path,
    # FIXME hack to get stdout logged
    logger,
  ):
    self.root_ctg            = root_ctg
    self.ctgfasta_path       = ctgfasta_path
    self.reads_ctg_bam_path  = reads_ctg_bam_path
    self.asmrootdir_path     = asmrootdir_path
    self.logger = logger

    self.tenxfq_paths = list(glob.glob(longranger_fqs_path + '/chnk*/files/*fastq'))

    self.fqdir_path = os.path.join(self.asmrootdir_path, 'fqs')
    util.mkdir_p(self.asmrootdir_path)
    util.mkdir_p(self.fqdir_path)

    self.ctg_size_map = util.get_fasta_sizes(self.ctgfasta_path)
 
  def assemble(self, local_asms, filt_ctgs=None):
    
    results = []
    for i, local_asm in enumerate(local_asms):
      if filt_ctgs and local_asm.link_ctg in filt_ctgs:
        self.logger.log('  - skipping filtered contig {}'.format(
          local_asm.link_ctg))
        continue
      contig_path = self._do_idba_assembly(local_asm)
      results.append((local_asm, contig_path))
    
    return results

  def _do_idba_assembly(self, local_asm):

    #continue
    self.logger.log('assembling with neighbor {}'.format(local_asm.link_ctg))
    self.logger.log('  - {} orig barcodes'.format(len(local_asm.bcode_set)))
    self.logger.log('  - {} downsampled barcodes'.format(len(local_asm.ds_bcode_set)))

    lrhintsfa_path = os.path.join(self.fqdir_path, 'lr-hints.{}.fa'.format(local_asm.uid))
    readsfa_path = os.path.join(self.fqdir_path, 'reads.{}.fa'.format(local_asm.uid))
    asmdir_path = os.path.join(self.asmrootdir_path, 'local-asm.{}'.format(local_asm.uid))

    # filter input reads
    get_bcode_reads(
      self.tenxfq_paths,
      readsfa_path,
      local_asm.ds_bcode_set,
    )
    # create long read hints
    hints = [(local_asm.root_ctg, local_asm.root_pos)]
    if local_asm.link_ctg:
      hints.append((local_asm.link_ctg, local_asm.link_pos))
    get_lrhints_fa(
      hints,
      self.ctgfasta_path,
      lrhintsfa_path,
    )

    #cmd = '{} -r {} -l {} -o {}'.format(
    cmd = '{} --maxk 60 -r {} -l {} -o {}'.format(
      idbabin_path,
      readsfa_path,
      lrhintsfa_path,
      asmdir_path,
    )
    pp = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
    retcode = pp.wait()
    print "::::", retcode
    contig_path = os.path.join(asmdir_path, 'contig.fa')
    if not os.path.exists(contig_path):
      if retcode != 0 and  "invalid insert distance" in pp.stderr.read():
        self.logger.log("idba_ud internal error; this is probably a bug in idba_ud, so we have to bail here")
        open(contig_path, "w")
      else:
        error_message = "assembly failed to produce contig.fa"
        self.logger.log(error_message)
        raise Exception()
    elif retcode != 0:
      self.logger.log(
        "something funky happened while running idba_ud (got error code "
        "{}) but there's not much we can do so continuing".format(retcode))
    # check contig.fa is valid
    try:
      fasta = pysam.FastaFile(contig_path)
    except Exception as e:
      self.logger.log('failed to produce correctly formatted contigs in {}'.format(contig_path))
      return None

    return contig_path

  def gen_local_cands(self):

    SMALL_CTG_SIZE = 1000
    COV_BIN_SIZE = 100
  
    fhandle = pysam.Samfile(self.reads_ctg_bam_path, 'rb')
    id_gen = util.IdGenerator()

    def get_barcode(read):
      filt_list = filter(lambda(k, v): k == 'BX', read.tags)
      if filt_list == []: 
        return None
      else:
        k, v = filt_list[0]
        return v
  
    def get_ctg_links(ctg, begin=None, end=None):
      bcodes = set()
      bcode_counts = Counter()
      link_ctg_counts = Counter()
      link_ctg_qnames = defaultdict(set)
      link_regions_map = defaultdict(lambda: ClusterTree(200, 1))
      link_reads_map = {}
      link_ctg_pos = {}
      cov_bin_counts = Counter()
      if begin == None or end == None:
        reads_iter = fhandle.fetch(ctg)
      else:
        assert None not in [begin, end]
        reads_iter = fhandle.fetch(ctg, begin, end)
      for i, read in enumerate(reads_iter):
        if read.is_unmapped:
          continue
        # skip low quality links for now
        if read.mapq < 10:
          continue
        bcode = get_barcode(read)
        if bcode == None:
          continue
        bcodes.add(bcode)
        bcode_counts[bcode] += 1
        tid = read.next_reference_id
        if tid != -1:
          p_ctg = fhandle.getrname(read.next_reference_id)    
          if p_ctg != ctg:
            link_ctg_counts[p_ctg] += 1
            link_ctg_qnames[p_ctg].add(read.qname)
            link_regions_map[p_ctg].insert(read.pos, read.aend, i)
            link_reads_map[i] = read
  
        # update coverage profile
        bin_idx = read.pos / COV_BIN_SIZE
        cov_bin_counts[bin_idx] += 1
        
      # filter for only links in which the paired end reads are clustered
      # near each other
      for p_ctg, regions in link_regions_map.items():
        regions = link_regions_map[p_ctg]
        (b, e, link_idxs) = max(regions.getregions(), key=lambda(x): len(x[2]))
        if len(link_idxs) < 3:
          del link_ctg_counts[p_ctg]
          del link_ctg_qnames[p_ctg]
        else:
          link_reads = map(lambda(i): link_reads_map[i], link_idxs)
          is_reverse = get_mode(map(lambda(r): r.is_reverse, link_reads))
          if is_reverse:
            pos = int(np.median(map(lambda(r): r.pos, link_reads)))
          else:
            pos = int(np.median(map(lambda(r): r.aend, link_reads)))
          link_ctg_pos[p_ctg] = (pos, is_reverse)

      return (
        bcodes,
        bcode_counts,
        link_ctg_counts,
        link_ctg_qnames,
        link_ctg_pos,
        cov_bin_counts,
      )

    def is_whack_coverage(cov_bin_counts):
      max_cov = max(cov_bin_counts.values())
      med_cov = np.median(cov_bin_counts.values())
      return max_cov > 20 * med_cov

    def get_mode(vs):
      return Counter(vs).most_common(1)[0][0]

    ( 
      root_bcode_set, 
      root_bcode_counts, 
      link_ctg_counts,
      link_ctg_qnames,
      root_link_ctg_pos,
      cov_bin_counts,
    ) = get_ctg_links(self.root_ctg)
  
    # don't locally assemble any small contigs with out of whack coverage
    if self.ctg_size_map[self.root_ctg] < SMALL_CTG_SIZE and is_whack_coverage(cov_bin_counts):
      max_cov = max(cov_bin_counts.values())
      med_cov = np.median(cov_bin_counts.values())
      self.logger.log('max coverage {} order of magnitude higher than median {}'.format(
        max_cov, med_cov))
      return []
  
    link_cand_ctgs = filter(
      lambda(c): (c != self.root_ctg and  link_ctg_counts[c] >= 3),
      link_ctg_counts,
    )

    self.logger.log('{} initial link candidates to examine'.format(
      len(link_cand_ctgs)))
    # examine recipricol links 
    link_ctgs = []
    # <ctg> : set(bcodes)
    link_ctg_bcodes = {}
    # <ctg> : Counter(bcodes)
    link_ctg_bcode_counts = {}
    # <ctg> : { <ctg> : pos }
    link_ctg_pos_map = {}
    for i, link_ctg in enumerate(link_cand_ctgs):
      #print 'link_ctg', link_ctg
      (
        link_ctg_bcodes[link_ctg],
        link_ctg_bcode_counts[link_ctg],
        _link_ctg_counts,
        _link_ctg_qnames,
        link_ctg_pos_map[link_ctg],
        _cov_bin_counts,
      ) = get_ctg_links(link_ctg)
      # link is not recripricol
      if _link_ctg_counts[self.root_ctg] < 3:
        continue
      # link contig has whacky coverage profile
      # FIXME consider removing this, may not be necessary
      if (
        self.ctg_size_map[link_ctg] < SMALL_CTG_SIZE and
        is_whack_coverage(_cov_bin_counts)
      ):
        continue
      # not enough barcode overlap
      if len(link_ctg_bcodes[link_ctg] & root_bcode_set) < 15:
        continue
      
      #print '  - added'
      link_ctgs.append(link_ctg)

    self.logger.log('  - {} pass reciprocal filtering'.format(
      len(link_ctgs)))

    local_asms = []
    for link_ctg in link_ctgs:
      if link_ctg == self.root_ctg:
        continue
      assert link_ctg_counts[link_ctg] >= 3
      link_bcode_counts = link_ctg_bcode_counts[link_ctg]
      i_bcode_set = link_ctg_bcodes[link_ctg] & root_bcode_set
      ds_bcode_set = i_bcode_set
      # downsample to barcodes with the most reads
      if len(i_bcode_set) > 400:
        ds_bcode_set = set(sorted(
          i_bcode_set,
          reverse=True,
          key=lambda(b): root_bcode_counts[b] + link_bcode_counts[b],
        )[:400])

      # must be recipricol link
      assert self.root_ctg in link_ctg_pos_map[link_ctg], \
        "recipricol link in link ctg {} to root {} not found".format(
          link_ctg,
          self.root_ctg,
        )
      local_asms.append(
        LocalAssembly(
          id_gen.get_next(),
          self.root_ctg,
          root_link_ctg_pos[link_ctg],
          link_ctg,
          link_ctg_pos_map[link_ctg][self.root_ctg],
          i_bcode_set,
          ds_bcode_set,
        )
      )

    def downsample(bcodes):
      if len(bcodes) > 400:
        return set(random.sample(root_bcode_set, 400))
      else:
        return set(bcodes)
  
    # do local reassembly of the root contig

    # if the contig is >20kb, create a local assembly at each end of the
    # root contig, otherwise just a single local assembly for the root
    root_size = self.ctg_size_map[self.root_ctg]
    if root_size > 2 * SEED_SELF_ASM_SIZE:
      self.logger.log('large root contig of size {}'.format(root_size))
      self.logger.log('  - generate head+tail local assemblies')
      # head
      (head_bcode_set, _, _, _, _, _,) = \
        get_ctg_links(self.root_ctg, 0, SEED_SELF_ASM_SIZE)
      begin_pos = (0, True)
      local_asms.append(
        LocalAssembly(
          id_gen.get_next(),
          self.root_ctg,
          begin_pos,
          None,
          None,
          head_bcode_set,
          downsample(head_bcode_set),
        )
      )
      # tail
      (tail_bcode_set, _, _, _, _, _,) = \
        get_ctg_links(self.root_ctg, root_size - SEED_SELF_ASM_SIZE, root_size)
      end_pos = (self.ctg_size_map[self.root_ctg], False)
      local_asms.append(
        LocalAssembly(
          id_gen.get_next(),
          self.root_ctg,
          end_pos,
          None,
          None,
          tail_bcode_set,
          downsample(tail_bcode_set),
        )
      )
    else:
      local_asms.append(
        LocalAssembly(
          id_gen.get_next(),
          self.root_ctg,
          None,
          None,
          None,
          root_bcode_set,
          downsample(root_bcode_set),
        )
      )
  
    return local_asms
  
def get_lrhints_fa(
  ctgs,
  infa_path,
  outfa_path,
):
  with open(outfa_path, 'w') as fout:
    ctg_fasta = pysam.FastaFile(infa_path)
    for ctg, pos_pre in ctgs:
      seq = str(ctg_fasta.fetch(ctg).upper())
      seq_sub = seq
      if pos_pre != None:
        pos, is_reverse = pos_pre
        if is_reverse:
          seq_sub = seq[pos+200:]
        else:
          seq_sub = seq[:max(0, pos-200)]
      # skip empty or uninformative long reads
      if len(seq_sub) < 200:
        continue
      for _ in xrange(10):
        fout.write('>{}\n{}\n'.format(ctg, seq_sub))

def get_bcode_reads(infq_paths, outfa_path, bcode_set):
  seen_set = set()
  with open(outfa_path, 'w') as fout:
    for fq_path in infq_paths:
      with FastqIndex(fq_path) as idx:
        for bcode in bcode_set & idx.bcode_set:
          for e in idx.get_reads(bcode):
            # tag qname with barcode
            nqname = '{}${}'.format(bcode, e.qname)
            fout.write('>{}\n'.format(nqname))
            fout.write('{}\n'.format(e.seq1))
            fout.write('>{}\n'.format(nqname))
            fout.write('{}\n'.format(e.seq2))
            seen_set.add(e.bcode)
  assert seen_set.issubset(bcode_set), 'not all barcodes loaded'
  return

