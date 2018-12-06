import os
import pysam
import subprocess
from collections import defaultdict, Counter
import glob
import random
import numpy as np
from bx.intervals.intersection import IntervalTree

from athena.mlib import util
from athena.mlib.fq_idx import FastqIndex

# NOTE must be in path
idbabin_path = 'idba_subasm'

SEED_SELF_ASM_SIZE = 10000
DS_SUBASM_COV = 100

SEED_SELF_ASM_SIZE = 4000
DS_SUBASM_COV = 60

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
    local_asm_cov,
  ):
    self.uid = uid
    self.root_ctg = root_ctg
    self.root_pos = root_pos
    self.link_ctg = link_ctg
    self.link_pos = link_pos
    self.bcode_set = bcode_set
    self.ds_bcode_set = ds_bcode_set
    self.local_asm_cov = local_asm_cov

  def __str__(self):
    return '{}_{}.{}'.format(self.root_ctg, self.link_ctg, self.uid)

class LocalAssembler(object):

  def __init__(
    self,
    root_ctg,
    ctgfasta_path,
    reads_ctg_bam_path,
    input_fqs,
    asmrootdir_path,
    # FIXME hack to get stdout logged
    logger,
  ):
    self.root_ctg            = root_ctg
    self.ctgfasta_path       = ctgfasta_path
    self.reads_ctg_bam_path  = reads_ctg_bam_path
    self.asmrootdir_path     = asmrootdir_path
    self.logger = logger

    self.tenxfq_paths = list(glob.glob(input_fqs))

    self.fqdir_path = os.path.join(self.asmrootdir_path, 'fqs')
    util.mkdir_p(self.asmrootdir_path)
    util.mkdir_p(self.fqdir_path)
    self.debugdir_path = os.path.join(self.asmrootdir_path, 'debug')
    util.mkdir_p(self.debugdir_path)

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
    self.logger.log('  - {}x estimated local coverage'.format(local_asm.local_asm_cov))

    lrhintsfa_path = os.path.join(self.fqdir_path, 'lr-hints.{}.fa'.format(local_asm.uid))
    seedsfa_path = os.path.join(self.fqdir_path, 'seeds.{}.fa'.format(local_asm.uid))
    readsfa_path = os.path.join(self.fqdir_path, 'reads.{}.ds.fa'.format(local_asm.uid))
    asmdir_path = os.path.join(self.asmrootdir_path, 'local-asm.{}'.format(local_asm.uid))
    fullreadsfa_path = os.path.join(self.fqdir_path, 'reads.{}.full.fa'.format(local_asm.uid))
    seed_info_path = os.path.join(self.asmrootdir_path, 'seeds-{}.txt'.format(local_asm.uid))
    
    with open(seed_info_path, 'w') as fout:
      fout.write(local_asm.root_ctg)
      if local_asm.root_pos != None:
        fout.write('\t'+str(local_asm.root_pos[0]))
      else:
        fout.write('\troot-asm')
      fout.write('\n')
      if local_asm.link_ctg == None:
        fout.write('no link-ctg\n')
      else:
        fout.write('{}\t{}\n'.format(local_asm.link_ctg, local_asm.link_pos[0]))

    # minimum required support based on estimated local_cov
    # NOTE do not allow anything to assemble with less than half of
    # min_support
    #min_support = max(2, int(local_asm.local_asm_cov / 10))
    # FIXME change back to /10
    min_support = max(2, int(local_asm.local_asm_cov / 20))
    self.logger.log('  - {} min_support required'.format(min_support))

    # filter input reads
    get_bcode_reads(
      self.tenxfq_paths,
      readsfa_path,
      local_asm.ds_bcode_set,
    )
    ## FIXME remove
    ## for now dump the non downsampled bcoded reads as well for debug
    #get_bcode_reads(
    #  self.tenxfq_paths,
    #  fullreadsfa_path,
    #  local_asm.bcode_set,
    #)
    # create long read hints
    hints = [(local_asm.root_ctg, local_asm.root_pos)]
    if local_asm.link_ctg:
      hints.append((local_asm.link_ctg, local_asm.link_pos))
    get_lrhints_fa(
      hints,
      self.ctgfasta_path,
      lrhintsfa_path,
      seedsfa_path,
      multiplicity=max(10, min_support+2),
    )

    cmd = '{} --num_threads 2 --min_support {} -r {} -l {} --seed_contig {} -o {}'.format(
      idbabin_path,
      min_support,
      readsfa_path,
      lrhintsfa_path,
      seedsfa_path,
      asmdir_path,
    )
    #cmd = '{} --maxk 80 -r {} -l {} -o {}'.format(
    #  idbabin_path,
    #  readsfa_path,
    #  lrhintsfa_path,
    #  asmdir_path,
    #)
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

    # compute read size
    # FIXME this is really hacky
    read_size = 0
    occ = 0
    for read in fhandle.fetch(self.root_ctg):
      if read.is_unmapped or read.is_secondary or read.is_supplementary:
        continue
      read_size = max(read_size, read.query_length)
      occ += 1
      if occ > 10:
        break

    def get_ctg_links(ctg, begin=None, end=None):
      bcodes = set()
      bcode_counts = Counter()
      link_ctg_counts = Counter()
      link_regions_map = defaultdict(lambda: IntervalTree(1, 1))
      link_reads_map = defaultdict(list)
      link_ctg_pos = {}
      cov_bin_counts = Counter()
      if begin == None or end == None:
        reads_iter = fhandle.fetch(ctg)
      else:
        assert None not in [begin, end]
        reads_iter = fhandle.fetch(ctg, begin, end)
      seen_set = set()
      for i, read in enumerate(reads_iter):
        if read.is_unmapped:
          continue
        # skip low quality links for now
        if read.mapq < 10:
          continue
        bcode = util.get_barcode(read)
        if bcode == None:
          continue
        optid = (read.pos, read.aend, bcode)
        # skip optical barcode duplicates
        if optid in seen_set:
          continue
        seen_set.add(optid)
        bcodes.add(bcode)
        bcode_counts[bcode] += 1
        tid = read.next_reference_id
        if tid != -1:
          p_ctg = fhandle.getrname(read.next_reference_id)    
          # save only links involving the root
          if p_ctg != ctg and self.root_ctg in [p_ctg, ctg]:
            link_ctg_counts[p_ctg] += 1
            link_regions_map[p_ctg].insert(read.pos, read.aend, read)
            link_reads_map[p_ctg].append(read)
  
        # update coverage profile
        bin_idx = read.pos / COV_BIN_SIZE
        cov_bin_counts[bin_idx] += 1
        
      # filter for only links in which the paired end reads are clustered
      # near each other
      for p_ctg, regions in link_regions_map.items():
        all_link_reads = link_reads_map[p_ctg]
        # find maximal link
        link_sets = map(
          lambda(r): regions.find(r.pos - 200, r.aend + 200),
          all_link_reads,
        )
        link_reads = max(link_sets, key=len)
        is_cand_chimera = len(all_link_reads) > 2 * len(link_reads)
        if len(link_reads) < 3 or is_cand_chimera:
          del link_ctg_counts[p_ctg]
        else:
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
      root_link_ctg_pos,
      cov_bin_counts,
    ) = get_ctg_links(self.root_ctg)
  
    # don't locally assemble any small contigs with out of whack coverage
    if self.ctg_size_map[self.root_ctg] < SMALL_CTG_SIZE and is_whack_coverage(cov_bin_counts):
      root_max_cov = max(cov_bin_counts.values())
      root_med_cov = np.median(cov_bin_counts.values())
      self.logger.log('root-ctg:{};filt-cov:True'.format(self.root_ctg))
      self.logger.log('max coverage {} order of magnitude higher than median {}'.format(
        root_max_cov, root_med_cov))
      return []
  
    link_cand_ctgs = filter(
      lambda(c): (c != self.root_ctg and  link_ctg_counts[c] >= 3),
      link_ctg_counts,
    )

    # truncate for a max of 200 checks
    num_orig_checks = len(link_cand_ctgs)
    truncate_checks = (num_orig_checks > 200)
    link_cand_ctgs = sorted(
      link_cand_ctgs,
      key=lambda(c): link_ctg_counts[c],
      reverse=True,
    )[:200]

    self.logger.log('{} initial link candidates to check'.format(
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
      (
        link_ctg_bcodes[link_ctg],
        link_ctg_bcode_counts[link_ctg],
        _link_ctg_counts,
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
      root_size = self.ctg_size_map[self.root_ctg]
      link_size = self.ctg_size_map[link_ctg]

      # refetch barcodes from locally linked region
      root_pos, is_rev = root_link_ctg_pos[link_ctg]
      if is_rev:
        rb, re = root_pos, min(root_pos + SEED_SELF_ASM_SIZE, root_size)
        _root_target_size = min(root_size - root_pos, SEED_SELF_ASM_SIZE)
      else:
        rb, re = max(0, root_pos - SEED_SELF_ASM_SIZE), root_pos
        _root_target_size = min(root_pos, SEED_SELF_ASM_SIZE)

      link_pos, is_rev = link_ctg_pos_map[link_ctg][self.root_ctg]
      if is_rev:
        lb, le = link_pos, min(link_pos + SEED_SELF_ASM_SIZE, link_size)
        _link_target_size = min(link_size - link_pos, SEED_SELF_ASM_SIZE)
      else:
        lb, le = max(0, link_pos - SEED_SELF_ASM_SIZE), link_pos
        _link_target_size = min(link_pos, SEED_SELF_ASM_SIZE)

      (_root_bcode_set, _root_bcode_counts, _, _, _,) = \
        get_ctg_links(self.root_ctg, rb, re)
      (_link_bcode_set, _link_bcode_counts, _, _, _,) = \
        get_ctg_links(link_ctg, lb, le)

      i_bcode_set = _link_bcode_set & _root_bcode_set
      ds_bcode_set, ds_num_reads = downsample_link_subassembly(
        _root_bcode_counts,
        _root_target_size,
        _link_bcode_counts,
        _link_target_size,
        read_size,
      )

      # estimate coverage of the barcodes for this local assembly
      local_asm_cov = 1. * ds_num_reads * read_size / (_root_target_size + _link_target_size)

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
          local_asm_cov,
        )
      )

    # truncate for a max of 50 local assemblies
    num_orig_asms = len(local_asms)
    truncate_asms = (num_orig_asms > 50)
    local_asms = sorted(
      local_asms,
      key=lambda(l): link_ctg_counts[l.link_ctg],
      reverse=True,
    )[:50]

    # do local reassembly of the root contig
    # if the contig is >20kb, create a local assembly at each end of the
    # root contig, otherwise just a single local assembly for the root
    root_size = self.ctg_size_map[self.root_ctg]
    if root_size > 1.5 * SEED_SELF_ASM_SIZE:
    #if root_size > 2 * SEED_SELF_ASM_SIZE:
      self.logger.log('large root contig of size {}'.format(root_size))
      self.logger.log('  - generate head+tail local assemblies')
      # head
      (head_bcode_set, head_root_bcode_counts, _, _, _,) = \
        get_ctg_links(self.root_ctg, 0, SEED_SELF_ASM_SIZE)
      ds_bcode_set, ds_num_reads = downsample_root_subassembly(
        head_root_bcode_counts,
        SEED_SELF_ASM_SIZE,
        read_size,
      )
      local_asm_cov = 1. * ds_num_reads * read_size / SEED_SELF_ASM_SIZE

      begin_pos = (0, True)
      local_asms.append(
        LocalAssembly(
          id_gen.get_next(),
          self.root_ctg,
          begin_pos,
          None,
          None,
          head_bcode_set,
          ds_bcode_set,
          local_asm_cov,
        )
      )
      # tail
      (tail_bcode_set, tail_root_bcode_counts, _, _, _,) = \
        get_ctg_links(self.root_ctg, root_size - SEED_SELF_ASM_SIZE, root_size)
      ds_bcode_set, ds_num_reads = downsample_root_subassembly(
        tail_root_bcode_counts,
        SEED_SELF_ASM_SIZE,
        read_size,
      )
      local_asm_cov = 1. * ds_num_reads * read_size / SEED_SELF_ASM_SIZE

      end_pos = (self.ctg_size_map[self.root_ctg], False)
      local_asms.append(
        LocalAssembly(
          id_gen.get_next(),
          self.root_ctg,
          end_pos,
          None,
          None,
          tail_bcode_set,
          ds_bcode_set,
          local_asm_cov,
        )
      )
    else:
      ds_bcode_set, num_reads = downsample_root_subassembly(
        root_bcode_counts,
        root_size,
        read_size,
      )
      local_asm_cov = 1. * num_reads * read_size / root_size
      local_asms.append(
        LocalAssembly(
          id_gen.get_next(),
          self.root_ctg,
          None,
          None,
          None,
          root_bcode_set,
          ds_bcode_set,
          local_asm_cov,
        )
      )
  
    fhandle.close()

    # debug logging message
    self.logger.log('root-ctg:{};numreads:{};checks:{};trunc-checks:{};asms:{};trunc-asms:{}'.format(
      self.root_ctg,
      sum(cov_bin_counts.values()),
      num_orig_checks,
      truncate_checks,
      num_orig_asms,
      truncate_asms,
    ))

    return local_asms
  
def get_lrhints_fa(
  ctgs,
  infa_path,
  outfa_path,
  seedsfa_path,
  multiplicity=10
):
  with open(outfa_path, 'w') as fout1, \
       open(seedsfa_path, 'w') as fout2:
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
      for _ in xrange(multiplicity):
        fout1.write('>{}\n{}\n'.format(ctg, seq_sub))
      fout2.write('>{}\n{}\n'.format(ctg, seq_sub))

def get_bcode_reads(infq_paths, outfa_path, bcode_set):
  seen_set = set()
  with open(outfa_path, 'w') as fout:
    for fq_path in infq_paths:
      with FastqIndex(fq_path) as idx:
        for bcode in bcode_set & idx.bcodes:
          for _, qname, lines in idx.get_reads(bcode):
            # tag qname with barcode
            nqname = '{}${}'.format(bcode, qname)
            seq = lines[1].strip()
            fout.write('>{}\n'.format(nqname))
            fout.write('{}\n'.format(seq))
          seen_set.add(bcode)
  assert seen_set.issubset(bcode_set), 'not all barcodes loaded'
  return

def downsample_link_subassembly(
  bcode_counts1,
  target_size1,
  bcode_counts2,
  target_size2,
  read_size,
):
  num_reads = 0
  i_bcodes = set(bcode_counts1.keys()) & set(bcode_counts2.keys())
  bcode_counts = Counter(
    {b: bcode_counts1[b] + bcode_counts2[b] for b in i_bcodes}
  )
  target_size = target_size1 + target_size2
  sub_bcodes = set()
  for bcode, _nr in bcode_counts.most_common():
    num_reads += _nr
    if 1. * num_reads * read_size / target_size < DS_SUBASM_COV:
      sub_bcodes.add(bcode)
    else:
      break
  return sub_bcodes, num_reads

def downsample_root_subassembly(bcode_counts, target_size, read_size):
  num_reads = 0
  sub_bcodes = set()
  for bcode, _nr in bcode_counts.most_common():
    num_reads += _nr
    if 1. * num_reads * read_size / target_size < DS_SUBASM_COV:
      sub_bcodes.add(bcode)
    else:
      break
  return sub_bcodes, num_reads

